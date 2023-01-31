# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adrian Mirza & Adarsh Kalikadien      #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #
import sys, os, glob

import morfeus
from morfeus import read_xyz, BiteAngle, ConeAngle, BuriedVolume, Dispersion, SASA
from morfeus.conformer import ConformerEnsemble
from morfeus.io import read_cclib, write_xyz
from morfeus.utils import convert_elements
import numpy as np
import pandas as pd

from molecular_graph import molecular_graph
from tools.utilities import dataframe_from_dictionary, calculate_distance
from dft_extraction import DFTExtractor


class Descriptors:
    """

    Class for calculating Morfeus and DFT descriptors. Works for Gaussian log files, CREST folders and xyz files.

    """
    def __init__(self, central_atom, path_to_workflow, output_type):
        self.central_atom = central_atom
        self.path_to_workflow = path_to_workflow
        self.supported_output_types = ['xyz', 'crest', 'gaussian']
        if output_type not in self.supported_output_types:
            raise ValueError(f'Output type {output_type} not supported. Please choose from {self.supported_output_types}.')
        self.output_type = output_type
        self.descriptor_df = None

    @staticmethod
    def _find_bidentate_ligand(elements, coordinates, geom_type):
        ligand_atoms, bidentate = molecular_graph(elements=elements, coords=coordinates, geom=geom_type)
        return ligand_atoms, bidentate

    def _merge_descriptor_dfs(self, old_descriptor_df, new_descriptor_df):
        old_descriptor_df = old_descriptor_df.merge(new_descriptor_df, on=['filename_tud',
                                                                            f"index_{self.central_atom}",
                                                                            "index_donor_max",
                                                                            "index_donor_min",
                                                                            f"element_{self.central_atom}",
                                                                            "element_donor_max",
                                                                            "element_donor_min"], how='left')
        return old_descriptor_df

    def set_output_type(self, new_output_type):
        """
        Set the output type of the descriptor calculator. This is used to determine how to read the files
        """
        if new_output_type not in self.supported_output_types:
            raise ValueError(f'Output type {new_output_type} not supported. Please choose from {self.supported_output_types}.')
        self.output_type = new_output_type

    def _calculate_steric_electronic_desc_morfeus(self, geom_type, solvent, dictionary, elements, coordinates):
        ligand_atoms, bidentate = self._find_bidentate_ligand(elements, coordinates, geom_type)
        # first index is the metal, second index is the bidentate ligand 1, third index is the bidentate ligand 2
        # morfeus indices start at 1, so add 1 to the indices
        metal_idx = bidentate[0] + 1
        bidentate_1_idx = bidentate[1] + 1
        bidentate_2_idx = bidentate[2] + 1

        # calculate steric descriptors
        dictionary["bite_angle"] = BiteAngle(coordinates, metal_idx, bidentate_1_idx,
                                             bidentate_2_idx).angle  # unit: degrees

        if geom_type == "BD" or geom_type == "SP":
            dictionary["cone_angle"] = ConeAngle(elements, coordinates, metal_idx).cone_angle
        else:
            try:
                a = [bidentate[0]]
                a.extend(ligand_atoms[bidentate[1]])
                a = list(np.sort(np.array(a)))

            except Exception:
                print('Molecular graph search failed, defaulting to manual search.')
                a = list(np.sort(np.array(bidentate)))

            diff = None
            for id, i in enumerate(a):
                if i == bidentate[0]:
                    diff = id
            elements_cone_angle = elements[a]
            coordinates_cone_angle = np.array(coordinates)[a]
            if diff is not None:
                dictionary["cone_angle"] = ConeAngle(elements_cone_angle, coordinates_cone_angle,
                                                     diff + 1).cone_angle
                # unit: degrees
            print('Cone angle calculation failed, defaulting to None. For complex:', complex)
            dictionary["cone_angle"] = None

        # determine max or min donor based on xTB charge of the donor atoms
        xtb_functional = 2  # indicate whether GFN1 or GFN2 is used for electronic descriptors
        if xtb_functional == 1:
            calculation_method = 'gfn1_xtb'
        else:
            calculation_method = 'gfn2_xtb'

        xtb = morfeus.XTB(elements, coordinates, solvent=solvent)
        atomic_charges = xtb.get_charges()  # determine min/max donor based on xTB charge of the donor atoms
        charge_bidentate_1 = atomic_charges[bidentate_1_idx]
        charge_bidentate_2 = atomic_charges[bidentate_2_idx]

        if charge_bidentate_1 > charge_bidentate_2:
            # this means that bidentate 1 is max donor (charge is less negative, so stronger donor)
            bidentate_max_donor_idx = bidentate_1_idx
            bidentate_min_donor_idx = bidentate_2_idx
        else:
            # this means that bidentate 2 is max donor
            bidentate_max_donor_idx = bidentate_2_idx
            bidentate_min_donor_idx = bidentate_1_idx

        # write indices and elements of metal center, max donor, and min donor to dictionary
        dictionary[f"index_{self.central_atom}"] = metal_idx
        dictionary["index_donor_max"] = bidentate_max_donor_idx
        dictionary["index_donor_min"] = bidentate_min_donor_idx
        element_symbols = convert_elements(elements, output='symbols')
        dictionary[f"element_{self.central_atom}"] = element_symbols[metal_idx - 1]
        dictionary["element_donor_max"] = element_symbols[bidentate_max_donor_idx - 1]
        dictionary["element_donor_min"] = element_symbols[bidentate_min_donor_idx - 1]

        # calculate distances between metal and donors, units: angstrom
        dictionary[f"distance_{self.central_atom}_max_donor_{self.output_type.lower}"] = calculate_distance(coordinates[metal_idx - 1], coordinates[bidentate_max_donor_idx - 1])
        dictionary[f"distance_{self.central_atom}_min_donor_{self.output_type.lower}"] = calculate_distance(coordinates[metal_idx - 1], coordinates[bidentate_min_donor_idx - 1])

        # calculate buried volume descriptors
        bv_metal_center = BuriedVolume(elements, coordinates, metal_idx, radius=3.5).fraction_buried_volume
        bv_max_donor = BuriedVolume(elements, coordinates, bidentate_1_idx, radius=3.5).fraction_buried_volume
        bv_min_donor = BuriedVolume(elements, coordinates, bidentate_2_idx, radius=3.5).fraction_buried_volume
        # unit: fraction of volume occupied by atoms within 3.5A of the metal
        dictionary[f"buried_volume_{self.central_atom}_3.5A"] = bv_metal_center
        dictionary["buried_volume_donor_max"] = bv_max_donor
        dictionary["buried_volume_donor_min"] = bv_min_donor

        buried_volume_for_quad_oct = BuriedVolume(elements, coordinates, metal_idx,
                                                  z_axis_atoms=bidentate_max_donor_idx,
                                                  xz_plane_atoms=[bidentate_min_donor_idx],
                                                  radius=3.5).octant_analysis()

        quadrants = buried_volume_for_quad_oct.quadrants['percent_buried_volume']
        octants = buried_volume_for_quad_oct.octants['percent_buried_volume']
        quadrant_dictionary = {1: 'NE', 2: 'NW', 3: 'SW', 4: 'SE'}
        octant_dictionary = {0: '+,+,+', 1: '-,+,+', 2: '-,-,+', 3: '+,-,+', 4: '+,-,-', 5: '-,-,-', 6: '-,+,-',
                             7: '+,+,-'}

        for quad_index in range(4):
            values = list(quadrants.values())
            dictionary[quadrant_dictionary[quad_index + 1] + "_quad"] = values[quad_index] / 100

        for oct_index in range(8):
            values = list(octants.values())
            dictionary[octant_dictionary[oct_index] + "_octant"] = values[oct_index] / 100

        bv_metal_4 = BuriedVolume(elements, coordinates, metal_idx, radius=4).fraction_buried_volume
        bv_metal_5 = BuriedVolume(elements, coordinates, metal_idx, radius=5).fraction_buried_volume
        bv_metal_6 = BuriedVolume(elements, coordinates, metal_idx, radius=6).fraction_buried_volume
        bv_metal_7 = BuriedVolume(elements, coordinates, metal_idx, radius=7).fraction_buried_volume

        dictionary[f"buried_volume_{self.central_atom}_4A"] = bv_metal_4
        dictionary[f"buried_volume_{self.central_atom}_5A"] = bv_metal_5
        dictionary[f"buried_volume_{self.central_atom}_6A"] = bv_metal_6
        dictionary[f"buried_volume_{self.central_atom}_7A"] = bv_metal_7

        # print(BuriedVolume(ce.elements, conformer.coordinates, bidentate[1]).print_report())
        # BuriedVolume(elements, coordinates, bidentate[1], radius=4).print_report()

        # calculate dispersion, SASA and electronic descriptors using morfeus
        # dispersion P_int unit: kcal^0.5 mol^-0.5
        dictionary[f"dispersion_p_int_{self.central_atom}_{calculation_method}"] = \
            Dispersion(elements, coordinates).atom_p_int[metal_idx]
        dictionary[f"dispersion_p_int_donor_max_{calculation_method}"] = \
            Dispersion(elements, coordinates).atom_p_int[bidentate_max_donor_idx]
        dictionary[f"dispersion_p_int_donor_min_{calculation_method}"] = \
            Dispersion(elements, coordinates).atom_p_int[bidentate_min_donor_idx]
        # SASA unit: A^2
        dictionary[f"sasa_{calculation_method}"] = SASA(elements, coordinates).area

        dictionary[f"ip_{calculation_method}"] = xtb.get_ip(corrected=True)  # unit: eV
        dipole = xtb.get_dipole()
        dictionary[f"dipole_{calculation_method}"] = np.sqrt(dipole.dot(dipole))  # unit: Debye dipole moment
        dictionary[f"ea_{calculation_method}"] = xtb.get_ea()  # unit: eV
        dictionary[f"electrofugality_{calculation_method}"] = xtb.get_global_descriptor(
            variety='electrofugality', corrected=True)  # unit: eV
        dictionary[f"nucleofugality_{calculation_method}"] = xtb.get_global_descriptor(
            variety='nucleofugality', corrected=True)  # unit: eV
        dictionary[f"nucleophilicity_{calculation_method}"] = xtb.get_global_descriptor(
            variety='nucleophilicity', corrected=True)  # unit: eV
        dictionary[f"electrophilicity_{calculation_method}"] = xtb.get_global_descriptor(
            variety='electrophilicity', corrected=True)  # unit: eV

        homo = xtb.get_homo()
        lumo = xtb.get_lumo()
        dictionary[f"HOMO_LUMO_gap_{calculation_method}"] = lumo - homo
        
        return dictionary, metal_idx, bidentate_max_donor_idx, bidentate_min_donor_idx

    def _calculate_dft_descriptors_from_log(self, log_file, metal_idx, bidentate_max_donor_idx, bidentate_min_donor_idx, dictionary):
        dft = DFTExtractor(log_file, metal_idx, bidentate_max_donor_idx, bidentate_min_donor_idx)

        # thermodynamic descriptors
        sum_electronic_and_free_energy, sum_electronic_and_enthalpy, zero_point_correction, entropy = dft.extract_thermodynamic_descriptors()
        dictionary[f"sum_electronic_and_free_energy_dft"] = sum_electronic_and_free_energy
        dictionary[f"sum_electronic_and_enthalpy_dft"] = sum_electronic_and_enthalpy
        dictionary[f"zero_point_correction_dft"] = zero_point_correction
        dictionary[f"entropy_dft"] = entropy

        # orbital occupations
        # donor with metal
        min_donor_metal_orbital_occupation, min_donor_metal_anti_orbital_occupation = dft.calculate_min_donor_metal_orbital_occupation(), dft.calculate_min_donor_metal_anti_orbital_occupation()
        dictionary[f"orbital_occupation_min_donor_{self.central_atom}_dft"] = min_donor_metal_orbital_occupation
        dictionary[f"anti_orbital_occupation_min_donor_{self.central_atom}_dft"] = min_donor_metal_anti_orbital_occupation
        max_donor_metal_orbital_occupation, max_donor_metal_anti_orbital_occupation = dft.calculate_max_donor_metal_orbital_occupation(), dft.calculate_max_donor_metal_anti_orbital_occupation()
        dictionary[f"orbital_occupation_max_donor_{self.central_atom}_dft"] = max_donor_metal_orbital_occupation
        dictionary[f"anti_orbital_occupation_max_donor_{self.central_atom}_dft"] = max_donor_metal_anti_orbital_occupation
        # donor with any other element
        min_donor_other_element_index_list, min_donor_other_orbital_occupation_list = dft.calculate_min_donor_other_orbital_occupation()
        if not min_donor_other_element_index_list is None and not min_donor_other_orbital_occupation_list is None:
            for i, (element_and_index, occupation) in enumerate(zip(min_donor_other_element_index_list, min_donor_other_orbital_occupation_list)):
                other_element = element_and_index[0]
                other_element_index = element_and_index[1]
                dictionary[f"orbital_occupation_min_donor_other_atom_{i + 1}_dft"] = occupation
                dictionary[f"orbital_occupation_min_donor_other_atom_{i + 1}_element_dft"] = other_element
                dictionary[f"orbital_occupation_min_donor_other_atom_{i + 1}_index_dft"] = other_element_index
        else:
            dictionary[f"orbital_occupation_min_donor_other_atom_1_dft"] = None
            dictionary[f"orbital_occupation_min_donor_other_atom_1_element_dft"] = None
            dictionary[f"orbital_occupation_min_donor_other_atom_1_index_dft"] = None

        min_donor_other_element_index_anti_bonding_list, min_donor_other_anti_orbital_occupation_list = dft.calculate_min_donor_other_anti_orbital_occupation()
        if not min_donor_other_element_index_anti_bonding_list is None and not min_donor_other_anti_orbital_occupation_list is None:
            for i, (element_and_index, occupation) in enumerate(zip(min_donor_other_element_index_anti_bonding_list, min_donor_other_anti_orbital_occupation_list)):
                other_element = element_and_index[0]
                other_element_index = element_and_index[1]
                dictionary[f"anti_orbital_occupation_min_donor_other_atom_{i + 1}_dft"] = occupation
                dictionary[f"anti_orbital_occupation_min_donor_other_atom_{i + 1}_element_dft"] = other_element
                dictionary[f"anti_orbital_occupation_min_donor_other_atom_{i + 1}_index_dft"] = other_element_index
        else:
            dictionary[f"anti_orbital_occupation_min_donor_other_atom_1_dft"] = None
            dictionary[f"anti_orbital_occupation_min_donor_other_atom_1_element_dft"] = None
            dictionary[f"anti_orbital_occupation_min_donor_other_atom_1_index_dft"] = None

        max_donor_other_element_index_list, max_donor_other_orbital_occupation_list = dft.calculate_max_donor_other_orbital_occupation()
        if not max_donor_other_element_index_list is None and not max_donor_other_orbital_occupation_list is None:
            for i, (element_and_index, occupation) in enumerate(zip(max_donor_other_element_index_list, max_donor_other_orbital_occupation_list)):
                other_element = element_and_index[0]
                other_element_index = element_and_index[1]
                dictionary[f"orbital_occupation_max_donor_other_atom_{i + 1}_dft"] = occupation
                dictionary[f"orbital_occupation_max_donor_other_atom_{i + 1}_element_dft"] = other_element
                dictionary[f"orbital_occupation_max_donor_other_atom_{i + 1}_index_dft"] = other_element_index
        else:
            dictionary[f"orbital_occupation_max_donor_other_atom_1_dft"] = None
            dictionary[f"orbital_occupation_max_donor_other_atom_1_element_dft"] = None
            dictionary[f"orbital_occupation_max_donor_other_atom_1_index_dft"] = None

        max_donor_other_element_index_anti_bonding_list, max_donor_other_anti_orbital_occupation_list = dft.calculate_max_donor_other_anti_orbital_occupation()
        if not max_donor_other_element_index_anti_bonding_list is None and not max_donor_other_anti_orbital_occupation_list is None:
            for i, (element_and_index, occupation) in enumerate(zip(max_donor_other_element_index_anti_bonding_list, max_donor_other_anti_orbital_occupation_list)):
                other_element = element_and_index[0]
                other_element_index = element_and_index[1]
                dictionary[f"anti_orbital_occupation_max_donor_other_atom_{i + 1}_dft"] = occupation
                dictionary[f"anti_orbital_occupation_max_donor_other_atom_{i + 1}_element_dft"] = other_element
                dictionary[f"anti_orbital_occupation_max_donor_other_atom_{i + 1}_index_dft"] = other_element_index
        else:
            dictionary[f"anti_orbital_occupation_max_donor_other_atom_1_dft"] = None
            dictionary[f"anti_orbital_occupation_max_donor_other_atom_1_element_dft"] = None
            dictionary[f"anti_orbital_occupation_max_donor_other_atom_1_index_dft"] = None

        # dipole moment
        dipole_moment = dft.calculate_dipole_moment()
        dictionary[f"dipole_moment_dft"] = dipole_moment

        # lone pair occupancy
        lone_pair_occupancy_min_donor, lone_pair_occupancy_max_donor = dft.calculate_donor_lone_pair_occupancy()
        dictionary["lone_pair_occupancy_min_donor_dft"] = lone_pair_occupancy_min_donor
        dictionary["lone_pair_occupancy_max_donor_dft"] = lone_pair_occupancy_max_donor

        # dispersion energy
        dispersion_energy = dft.calculate_dispersion_energy()
        dictionary["dispersion_energy_dft"] = dispersion_energy

        # NBO charges
        metal_nbo_charge, min_donor_nbo_charge, max_donor_nbo_charge = dft.calculate_natural_charges()
        dictionary[f"nbo_charge_{self.central_atom}_dft"] = metal_nbo_charge
        dictionary[f"nbo_charge_min_donor_dft"] = min_donor_nbo_charge
        dictionary[f"nbo_charge_max_donor_dft"] = max_donor_nbo_charge

        # mulliken charges
        metal_mulliken_charge, min_donor_mulliken_charge, max_donor_mulliken_charge = dft.calculate_mulliken_charges()
        dictionary[f"mulliken_charge_{self.central_atom}_dft"] = metal_mulliken_charge
        dictionary[f"mulliken_charge_min_donor_dft"] = min_donor_mulliken_charge
        dictionary[f"mulliken_charge_max_donor_dft"] = max_donor_mulliken_charge

        # other electronic descriptors
        homo_energy, lumo_energy, homo_lumo_gap, hardness, softness, electronegativity, electrophilicity = dft.calculate_electronic_descriptors()
        dictionary["homo_energy_dft"] = homo_energy
        dictionary["lumo_energy_dft"] = lumo_energy
        dictionary["homo_lumo_gap_dft"] = homo_lumo_gap
        dictionary["hardness_dft"] = hardness
        dictionary["softness_dft"] = softness
        dictionary["electronegativity_dft"] = electronegativity
        dictionary["electrophilicity_dft"] = electrophilicity

        return dictionary

    def calculate_morfeus_descriptors(self, geom_type, solvent=None, printout=False):
        if self.output_type.lower() == 'xyz':
            complexes_to_calc_descriptors = glob.glob(os.path.join(self.path_to_workflow, '*.xyz'))
            dictionary_for_properties = {}

            # try:
            for metal_ligand_complex in complexes_to_calc_descriptors:
                properties = {}

                base_with_extension = os.path.basename(metal_ligand_complex)
                split_base = os.path.splitext(base_with_extension)
                filename = split_base[0]
                properties['filename_tud'] = filename

                elements, coordinates = read_xyz(metal_ligand_complex)
                properties, metal_idx, bidentate_max_donor_idx, bidentate_min_donor_idx = self._calculate_steric_electronic_desc_morfeus(geom_type=geom_type, solvent=solvent, dictionary=properties, elements=elements, coordinates=coordinates)
                dictionary_for_properties[os.path.basename(os.path.normpath(metal_ligand_complex[:-4]))] = properties

            new_descriptor_df = dataframe_from_dictionary(dictionary_for_properties)

            if printout:
                print(new_descriptor_df.to_markdown())

            # merge descriptor dataframes if they already exist
            if self.descriptor_df is None:
                self.descriptor_df = new_descriptor_df
            else:
                self.descriptor_df = self._merge_descriptor_dfs(self.descriptor_df, new_descriptor_df)

        elif self.output_type.lower() == 'crest':
            complexes_to_calc_descriptors = glob.glob(os.path.join(self.path_to_workflow, 'CREST', '*'))
            dictionary_for_conformer_properties = {}
            for complex in complexes_to_calc_descriptors:
                conformer_properties = {}
                ce = None
                try:
                    ce = ConformerEnsemble.from_crest(complex)
                except Exception as e:
                    print("Descriptor calculation failed for this complex:", complex)
                    print(e)
                    conformer_properties['filename_tud'] = os.path.basename(os.path.normpath(complex))
                    continue

                if ce is not None:
                    # ce.generate_mol()
                    ce.prune_energy()
                    ce.sort()
                    for conformer in ce:
                        elements, coordinates = ce.elements, conformer.coordinates
                        conformer.properties, metal_idx, bidentate_max_donor_idx, bidentate_min_donor_idx = self._calculate_steric_electronic_desc_morfeus(geom_type=geom_type, solvent=solvent, dictionary=conformer.properties, elements=elements, coordinates=coordinates)
                    # all descriptors calculated, now we can write the filaname and boltzman statistics to the dictionary
                    conformer_properties['filename_tud'] = os.path.basename(os.path.normpath(complex))

                    columns_to_exclude = [f"index_{self.central_atom}", "index_donor_max", "index_donor_min", f"element_{self.central_atom}", "element_donor_max", "element_donor_min"]
                    for key in [k for k in ce.get_properties().keys() if k in columns_to_exclude]:
                        # check if indexing property is the same across all conformers, then select property of first
                        # conformer
                        for property in ce.get_properties()[key]:
                            if property != ce.get_properties()[key][0]:
                                print(f"BE AWARE: Indexing property {key} is not the same across all conformers. Taking property "
                                      f"of first conformer for {os.path.basename(os.path.normpath(complex))}.")
                        conformer_properties[key] = ce.get_properties()[key][0]

                    # boltzmann averaging
                    for key in [k for k in ce.get_properties().keys() if k not in columns_to_exclude]:
                        conformer_properties[f"{key}_boltzmann_average"] = ce.boltzmann_statistic(key)
                        conformer_properties[f"{key}_boltzmann_std"] = ce.boltzmann_statistic(key, statistic='std')
                        conformer_properties[f"{key}_boltzmann_variance"] = ce.boltzmann_statistic(key, statistic='var')
                        conformer_properties[f"{key}_Emin_conformer"] = ce.get_properties()[key][0]
                dictionary_for_conformer_properties[os.path.basename(os.path.normpath(complex))] = conformer_properties

            new_descriptor_df = dataframe_from_dictionary(dictionary_for_conformer_properties)
            if printout:
                print(new_descriptor_df.to_markdown())

            if self.descriptor_df is None:
                self.descriptor_df = new_descriptor_df
            else:
                self.descriptor_df = self._merge_descriptor_dfs(self.descriptor_df, new_descriptor_df)

        else:
            raise ValueError(f'Output type {self.output_type} not supported. Please choose from {self.supported_output_types}.')

    def calculate_dft_descriptors_from_log(self, geom_type, solvent=None, extract_xyz_from_log=False, printout=False):
        # get all log files
        complexes_to_calc_descriptors = glob.glob(os.path.join(self.path_to_workflow, '*.log'))
        dictionary_for_properties = {}

        # first calculate morfeus descriptors in same way as for xyz files using cclib
        for metal_ligand_complex in complexes_to_calc_descriptors:
            properties = {}
            base_with_extension = os.path.basename(metal_ligand_complex)
            split_base = os.path.splitext(base_with_extension)
            filename = split_base[0]
            properties['filename_tud'] = filename

            elements, coordinates = read_cclib(metal_ligand_complex)
            coordinates = coordinates[-1]  # morfeus descriptors are calculated for last geometry in log file
            elements = np.array(elements)
            properties, metal_idx, bidentate_max_donor_idx, bidentate_min_donor_idx = self._calculate_steric_electronic_desc_morfeus(geom_type=geom_type, solvent=solvent, dictionary=properties, elements=elements, coordinates=coordinates)

            # calculate DFT descriptors from Gaussian log file
            # get indices of bidentate ligands and metal for descriptor calculation class
            properties = self._calculate_dft_descriptors_from_log(metal_ligand_complex, metal_idx, bidentate_max_donor_idx, bidentate_min_donor_idx, properties)

            # write xyz for log file
            if extract_xyz_from_log:
                xyz_filename = metal_ligand_complex[:-4] + '_DFT.xyz'
                write_xyz(os.path.join(self.path_to_workflow, xyz_filename), elements, coordinates)

            # for property in properties.keys():
            dictionary_for_properties[os.path.basename(os.path.normpath(metal_ligand_complex[:-4]))] = properties

        new_descriptor_df = dataframe_from_dictionary(dictionary_for_properties)

        if printout:
            print(new_descriptor_df.to_markdown())

        if self.descriptor_df is None:
            self.descriptor_df = new_descriptor_df
        else:
            self.descriptor_df = self._merge_descriptor_dfs(self.descriptor_df, new_descriptor_df)


if __name__ == "__main__":
    # descriptors = Descriptors(central_atom='Rh', path_to_workflow=os.path.join(os.getcwd(), 'Workflow'), output_type='xyz')
    # descriptors.calculate_morfeus_descriptors(geom_type='BD')
    # descriptors.descriptor_df.to_csv('descriptors.csv', index=False)

    conformer_descriptors = Descriptors(central_atom='Rh', path_to_workflow=os.path.join(os.getcwd(), 'Workflow'), output_type='gaussian')
    conformer_descriptors.calculate_dft_descriptors_from_log(geom_type='BD', solvent=None, extract_xyz_from_log=True, printout=False)
    conformer_descriptors.descriptor_df.to_csv('DFT_descriptors.csv', index=False)
    # conformer_descriptors.calculate_morfeus_descriptors(geom_type='BD')
    # # conformer_descriptors.descriptor_df.to_csv('conformer_descriptors.csv', index=False)
    # conformer_descriptors.set_output_type('xyz')
    # conformer_descriptors.calculate_morfeus_descriptors(geom_type='BD')
    # conformer_descriptors.descriptor_df.to_csv('conformer_descriptors', index=False)

