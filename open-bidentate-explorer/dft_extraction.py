# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #
import os
import glob
import re
import numpy as np
from tools.utilities import get_bonded_atoms
from morfeus.io import read_cclib, get_xyz_string, read_xyz
from morfeus.utils import convert_elements
import cclib
from cclib.parser import ccopen


class DFTExtractor(object):
    def __init__(self, log_file, metal_center_idx, min_donor_idx, max_donor_idx, metal_adduct='pristine'):
        self.log_file = log_file

        supported_metal_adducts = ['pristine', 'acetonitrile', 'nbd']  # norbornadiene is placed at bottom of xyz file, so it is a useful pointer for quadrant analysis
        if metal_adduct.lower() not in supported_metal_adducts:
            raise ValueError(
                f"Metal adduct {metal_adduct} not supported. Please choose from {supported_metal_adducts}.")
        self.metal_adduct = metal_adduct.lower()

        # use cclib parser
        self.parser = ccopen(log_file)
        self.data = self.parser.parse()
        self.meta_data = self.data.metadata
        if not self.check_normal_termination():
            print(f"Warning: {log_file} did not terminate normally!")

        self.atom_charges_dict = self.data.atomcharges
        # use morfeus parser (in case it's needed)
        self.elements, self.coordinates = read_cclib(log_file)
        if not len(self.coordinates[-1]) == 3:  # if this is true, there is only 1 coordinates array
            self.coordinates = self.coordinates[-1]  # else morfeus descriptors are calculated for last geometry in log file
        self.elements = np.array(self.elements)

        # get index of nbd C in back of molecule which we use as indicator for z-axis
        if self.metal_adduct == 'nbd':
            self.carbon_back_nbd_idx = np.where(self.coordinates == self.coordinates[-15])[0][0]
            # we need the xyz string to get the hydrogens bonded to this carbon
            self.xyz_string = ""
            coordinates = np.array(self.coordinates).reshape(-1, len(self.elements), 3)
            for coord in coordinates:
                self.xyz_string = get_xyz_string(convert_elements(self.elements, 'symbols'), coord)

        # idx's come from morfeus, so they start at 1, subtract 1 to get the correct index for cclib
        self.metal_center_idx = metal_center_idx - 1
        self.min_donor_idx = min_donor_idx - 1
        self.max_donor_idx = max_donor_idx - 1
        self.elements = convert_elements(self.elements, output='symbols')
        self.metal_center_element = self.elements[self.metal_center_idx]
        self.min_donor_element = self.elements[self.min_donor_idx]
        self.max_donor_element = self.elements[self.max_donor_idx]

        # read the log file for extraction
        with open(log_file) as f:
            self.log_file_text = f.readlines()

    def check_normal_termination(self):
        success = self.meta_data['success']
        return success

    def check_nbd_back_carbon(self):
        # ToDo: build some more checks in here before returning the idx
        if self.metal_adduct == 'nbd':
            if len(self.get_hydrogens_bonded_to_carbon_back_nbd()) == 2:
                return self.carbon_back_nbd_idx
        return None

    def get_hydrogens_bonded_to_carbon_back_nbd(self):
        if self.metal_adduct == 'nbd':
            # get the hydrogen indices
            hydrogen_indices = get_bonded_atoms(self.xyz_string, int(self.carbon_back_nbd_idx), 1)  # 1 is the hydrogen atom
            # this gives a list of atomic number and index, we only want the index
            if len(hydrogen_indices) == 2:
                return [hydrogen_indices[0][1] - 1, hydrogen_indices[1][1] - 1]
            else:
                print(f"Warning: {self.log_file} does not have 2 hydrogens bonded to the carbon in the back of the nbd molecule. This is needed for dihedral angle calculation. Skipping this molecule.")
                return None
        else:
            return None

    def extract_time(self):
        wall_time = self.meta_data['wall_time']
        cpu_time = self.meta_data['cpu_time']
        return wall_time, cpu_time

    def extract_homo_lumo_gap(self):
        homo = self.data.homos[0]
        energies = self.data.moenergies[0]
        lumo_energy = energies[homo+1]
        homo_energy = energies[homo]
        homo_lumo_gap = lumo_energy - homo_energy
        return homo_energy, lumo_energy, homo_lumo_gap

    def extract_natural_charges(self):
        # NBO charge
        natural_charges = self.atom_charges_dict['natural']
        return natural_charges

    def extract_mulliken_charges(self):
        atomic_charges = self.atom_charges_dict['mulliken']
        return atomic_charges

    def extract_dipole_moment(self):
        dipole_moment = self.data.moments
        return dipole_moment

    def extract_dispersion_energies(self):
        dispersion_energies = self.data.dispersionenergies
        return dispersion_energies

    # def extract_natural_orbital_coefficients(self):
    #     natural_orbital_coefficients = self.data.nocoeffs
    #     return natural_orbital_coefficients

    def extract_donor_metal_antibonding_occupation(self, donor_element, donor_idx, metal_center_element, metal_center_idx):
        data = self.log_file_text
        donor_lines = []
        nbo_indices = []
        for line_index, line in enumerate(data):
            if "NBO                        Occupancy    Energy   (geminal,vicinal,remote)" in line:
                nbo_indices.append(line_index)
        end_nbo_search_indices = []
        for line_index, line in enumerate(data):
            # every nbo summary ends with two lines of "       -------------------------------"
            if re.search("[ ]{7}-------------------------------", line):
                end_nbo_search_indices.append(line_index)
        # from each pair we only need the first one to mark the end of the nbo summary
        end_nbo_search_indices = [end_nbo_search_indices[i] for i in range(0, len(end_nbo_search_indices), 2)]

        # antibonding orbital for min donor
        for nr_of_nbo_summary_found, nbo_index in enumerate(nbo_indices):
            search_this_area = data[nbo_index:end_nbo_search_indices[
                nr_of_nbo_summary_found + 1]] if nr_of_nbo_summary_found < len(nbo_indices) - 1 else data[nbo_index:]
            for line in search_this_area:
                pattern = "BD\*\(\s+1\)" + f" {donor_element}\s+" + f"{donor_idx + 1}\s" + "[-]" + f"{metal_center_element}\s+" + f"{metal_center_idx + 1}" + ".*[\d\s!@#\\$%\\^\\&*\\)\\(+=._-]$"
                if re.search(pattern, line):
                    donor_lines.append(line)
                    # print(line)
                second_pattern = "BD\*\(\s+1\)" + f"{metal_center_element}\s+" + f"{metal_center_idx + 1}\s" + "[-]" + f"\s+{donor_element}\s+" + f"{donor_idx + 1}" + ".*[\d\s!@#\\$%\\^\\&*\\)\\(+=._-]$"
                if re.search(second_pattern, line):
                    donor_lines.append(line)
                    # print(line)

        # only digits are important
        donor_lines = [re.sub('[^\d\.]', ' ', line_).strip().split() for line_ in donor_lines]

        if len(donor_lines) == 0:
            return None
        else:
            return donor_lines[-1][4]

    def extract_donor_metal_bonding_occupation(self, donor_element, donor_idx, metal_center_element, metal_center_idx):
        data = self.log_file_text
        donor_lines = []
        nbo_indices = []
        for line_index, line in enumerate(data):
            if "NBO                        Occupancy    Energy   (geminal,vicinal,remote)" in line:
                nbo_indices.append(line_index)
        end_nbo_search_indices = []
        for line_index, line in enumerate(data):
            # every nbo summary ends with two lines of "       -------------------------------"
            if re.search("[ ]{7}-------------------------------", line):
                end_nbo_search_indices.append(line_index)
            # from each pair we only need the first one to mark the end of the nbo summary
        end_nbo_search_indices = [end_nbo_search_indices[i] for i in range(0, len(end_nbo_search_indices), 2)]

        # bonding orbital for min donor
        for nr_of_nbo_summary_found, nbo_index in enumerate(nbo_indices):
            search_this_area = data[nbo_index:end_nbo_search_indices[
                nr_of_nbo_summary_found]] if nr_of_nbo_summary_found < len(nbo_indices) else data[nbo_index:]
            for line in search_this_area:
                pattern = "BD\s\(\s+1\)" + f" {donor_element}\s+" + f"{donor_idx + 1}\s" + "[-]" + f"{metal_center_element}\s+" + f"{metal_center_idx + 1}" + ".*[\d\s!@#\\$%\\^\\&*\\)\\(+=._-]$"
                if re.search(pattern, line):
                    donor_lines.append(line)
                second_pattern = "BD\s\(\s+1\)" + f"{metal_center_element}\s+" + f"{metal_center_idx + 1}\s" + "[-]" + f"\s+{donor_element}\s+" + f"{donor_idx + 1}" + ".*[\d\s!@#\\$%\\^\\&*\\)\\(+=._-]$"
                if re.search(second_pattern, line):
                    donor_lines.append(line)

        # only digits are important
        donor_lines = [re.sub('[^\d\.]', ' ', line_).strip().split() for line_ in donor_lines]

        # we are only interested in the last iteration, so we take the last element of the list
        if len(donor_lines) == 0:
            return None
        else:
            return donor_lines[-1][4]

    def extract_donor_other_element_antibonding_occupation(self, donor_element, donor_idx):
        data = self.log_file_text
        donor_lines = []
        nbo_indices = []
        for line_index, line in enumerate(data):
            if "NBO                        Occupancy    Energy   (geminal,vicinal,remote)" in line:
                nbo_indices.append(line_index)
        end_nbo_search_indices = []
        for line_index, line in enumerate(data):
            # every nbo summary ends with two lines of "       -------------------------------"
            if re.search("[ ]{7}-------------------------------", line):
                end_nbo_search_indices.append(line_index)
            # from each pair we only need the first one to mark the end of the nbo summary
        end_nbo_search_indices = [end_nbo_search_indices[i] for i in range(0, len(end_nbo_search_indices), 2)]
        # bonding orbital for min donor
        for nr_of_nbo_summary_found, nbo_index in enumerate(nbo_indices):
            search_this_area = data[nbo_index:end_nbo_search_indices[
                nr_of_nbo_summary_found]] if nr_of_nbo_summary_found < len(end_nbo_search_indices) else data[nbo_index:]

            for line in search_this_area:
                pattern = "BD\*\(\s+1\)" + f" {donor_element}\s+" + f"{donor_idx + 1}\s" + "[-]" + "\s[CONHP]\s+" + ".*[\d\s!@#\\$%\\^\\&*\\)\\(+=._-]$"
                if re.search(pattern, line):
                    donor_lines.append(line)
                second_pattern = "BD\*\(\s+1\)" + " [CONHP]\s+" + "\d{1,2}\s" + "[-]" + f"\s+{donor_element}\s+" + f"{donor_idx + 1}" + ".*[\d\s!@#\\$%\\^\\&*\\)\\(+=._-]$"
                if re.search(second_pattern, line):
                    donor_lines.append(line)

        # if donor is P, there are 3 other bonds, so from the last 3 lines we take the values
        if donor_element == "P":
            donor_lines = donor_lines[-3:]
        elif donor_element == "N":
            donor_lines = donor_lines[-2:]
        else:
            raise ValueError("Donor element not supported")

        # get other element and index of that element
        pattern = "[CONHP]\s+\d{1,2}"
        other_element_and_index = [re.findall(pattern, line_) for line_ in donor_lines]
        flat_list = [item for sublist in other_element_and_index for item in sublist]
        # remove the donor element and idx from the list
        pattern = f"{donor_element}\s+{donor_idx + 1}"
        other_element_and_index = [i.split() for i in flat_list if not re.search(pattern, str(i))]

        # get occupation values by getting only digits from the lines
        occupation_values = [re.sub('[^\d\.]', ' ', line_).strip().split() for line_ in donor_lines]
        occupation_values = [value[4] for value in occupation_values]

        if len(donor_lines) != 0:
            if len(occupation_values) == len(other_element_and_index):
                # other element and index is a list of lists for each element, occupation_values is a list
                return other_element_and_index, occupation_values
            else:
                raise ValueError("Length of other element and index list and occupation value list is not equal")
        else:
            return None, None

    def extract_donor_other_element_bonding_occupation(self, donor_element, donor_idx):
        data = self.log_file_text
        donor_lines = []
        nbo_indices = []
        for line_index, line in enumerate(data):
            if "NBO                        Occupancy    Energy   (geminal,vicinal,remote)" in line:
                nbo_indices.append(line_index)
        end_nbo_search_indices = []
        for line_index, line in enumerate(data):
            # every nbo summary ends with two lines of "       -------------------------------"
            if re.search("[ ]{7}-------------------------------", line):
                end_nbo_search_indices.append(line_index)
            # from each pair we only need the first one to mark the end of the nbo summary
        end_nbo_search_indices = [end_nbo_search_indices[i] for i in range(0, len(end_nbo_search_indices), 2)]
        # bonding orbital for min donor
        for nr_of_nbo_summary_found, nbo_index in enumerate(nbo_indices):
            search_this_area = data[nbo_index:end_nbo_search_indices[
                nr_of_nbo_summary_found]] if nr_of_nbo_summary_found < len(end_nbo_search_indices) else data[nbo_index:]

            for line in search_this_area:
                pattern = "BD\s\(\s+1\)" + f" {donor_element}\s+" + f"{donor_idx + 1}\s" + "[-]" + "\s[CONHP]\s+" + ".*[\d\s!@#\\$%\\^\\&*\\)\\(+=._-]$"
                if re.search(pattern, line):
                    donor_lines.append(line)
                second_pattern = "BD\s\(\s+1\)" + " [CONHP]\s+" + "\d{1,2}\s" + "[-]" + f"\s+{donor_element}\s+" + f"{donor_idx + 1}" + ".*[\d\s!@#\\$%\\^\\&*\\)\\(+=._-]$"
                if re.search(second_pattern, line):
                    donor_lines.append(line)

        # if donor is P, there are 3 other bonds, so from the last 3 lines we take the values
        if donor_element == "P":
            donor_lines = donor_lines[-3:]
        elif donor_element == "N":
            donor_lines = donor_lines[-2:]
        else:
            raise ValueError("Donor element not supported")

        # get other element and index of that element
        pattern = "[CONHP]\s+\d{1,2}"
        other_element_and_index = [re.findall(pattern, line_) for line_ in donor_lines]
        flat_list = [item for sublist in other_element_and_index for item in sublist]
        # remove the donor element and idx from the list
        pattern = f"{donor_element}\s+{donor_idx + 1}"
        other_element_and_index = [i.split() for i in flat_list if not re.search(pattern, str(i))]

        # get occupation values by getting only digits from the lines
        occupation_values = [re.sub('[^\d\.]', ' ', line_).strip().split() for line_ in donor_lines]
        occupation_values = [value[4] for value in occupation_values]

        if len(donor_lines) != 0:
            if len(occupation_values) == len(other_element_and_index):
                # other element and index is a list of lists for each element, occupation_values is a list
                return other_element_and_index, occupation_values
            else:
                raise ValueError("Length of other element and index list and occupation value list is not equal")
        else:
            return None, None

    def extract_donor_lone_pair_occupancy(self, atom_type_min, atom_type_max):
        # Count nr of times " NATURAL POPULATIONS:  Natural atomic orbital occupancies" appears and take the last one
        count = 0
        count2 = 0
        # with open(self.log_file) as file:  # Use file to refer to the file object
        data = self.log_file_text
        for line_index, line in enumerate(data):
            if "NATURAL POPULATIONS:  Natural atomic orbital occupancies" in line:
                count += 1
                if count == 3:
                    counting_line = line_index + 4
            if "electrons found in the effective core potential" in line:
                count2 += 1
                if count2 == 3:
                    counting_line_ = line_index - 1

                # if
                # final_line =
        # file.close()
        occupancy_raw_data = data[counting_line: counting_line_]

        occupancy_final_data = []
        orbital_dictionary = {}

        for occupancy in occupancy_raw_data:
            occupancy = occupancy.strip("\n")
            occupancy = occupancy.strip("\r")
            occupancy = occupancy.split()
            if occupancy != []:
                occupancy_final_data.append((occupancy))

        for orbital in occupancy_final_data:
            if int(orbital[2]) not in orbital_dictionary.keys():
                orbital_dictionary[int(orbital[2])] = {orbital[3] + '_' + orbital[4] + orbital[5]: float(orbital[6])}
            else:
                orbital_dictionary[int(orbital[2])][orbital[3] + '_' + orbital[4] + orbital[5]] = float(orbital[6])
        # divided by 2 since 2 electrons form the lone pair and the 3s orbital

        returns = []
        if atom_type_min == 'N':
            returns.append(orbital_dictionary[self.min_donor_idx + 1]['S_Val(2S)'])

        else:
            returns.append(orbital_dictionary[self.min_donor_idx + 1]['S_Val(3S)'])

        if atom_type_max == 'N':
            returns.append(orbital_dictionary[self.max_donor_idx + 1]['S_Val(2S)'])

        else:
            returns.append(orbital_dictionary[self.max_donor_idx + 1]['S_Val(3S)'])

        return returns[0], returns[1]

    def extract_thermodynamic_descriptors(self):
        # extract energy, enthalpy, entropy, free energy in hartree
        sum_electronic_and_free_energy = self.data.freeenergy
        sum_electronic_and_enthalpy = self.data.enthalpy
        zero_point_correction = self.data.zpve
        entropy = self.data.entropy
        return sum_electronic_and_free_energy, sum_electronic_and_enthalpy, zero_point_correction, entropy

    def calculate_min_donor_metal_orbital_occupation(self):
        metal_min_donor_occupancy = self.extract_donor_metal_bonding_occupation(self.min_donor_element, self.min_donor_idx, self.metal_center_element, self.metal_center_idx)
        return metal_min_donor_occupancy

    def calculate_min_donor_metal_anti_orbital_occupation(self):
        metal_min_donor_anti_occupancy = self.extract_donor_metal_antibonding_occupation(self.min_donor_element, self.min_donor_idx, self.metal_center_element, self.metal_center_idx)
        return metal_min_donor_anti_occupancy

    def calculate_max_donor_metal_orbital_occupation(self):
        metal_max_donor_occupancy = self.extract_donor_metal_bonding_occupation(self.max_donor_element, self.max_donor_idx, self.metal_center_element, self.metal_center_idx)
        return metal_max_donor_occupancy

    def calculate_max_donor_metal_anti_orbital_occupation(self):
        metal_max_donor_anti_occupancy = self.extract_donor_metal_antibonding_occupation(self.max_donor_element, self.max_donor_idx, self.metal_center_element, self.metal_center_idx)
        return metal_max_donor_anti_occupancy

    def calculate_min_donor_other_orbital_occupation(self):
        other_element_index_list, other_element_occupancy = self.extract_donor_other_element_bonding_occupation(self.min_donor_element, self.min_donor_idx)
        return other_element_index_list, other_element_occupancy

    def calculate_min_donor_other_anti_orbital_occupation(self):
        other_element_index_list, other_element_occupancy = self.extract_donor_other_element_antibonding_occupation(self.min_donor_element, self.min_donor_idx)
        return other_element_index_list, other_element_occupancy

    def calculate_max_donor_other_orbital_occupation(self):
        other_element_index_list, other_element_occupancy = self.extract_donor_other_element_bonding_occupation(self.max_donor_element, self.max_donor_idx)
        return other_element_index_list, other_element_occupancy

    def calculate_max_donor_other_anti_orbital_occupation(self):
        other_element_index_list, other_element_occupancy = self.extract_donor_other_element_antibonding_occupation(self.max_donor_element, self.max_donor_idx)
        return other_element_index_list, other_element_occupancy

    def calculate_donor_lone_pair_occupancy(self):
        min_donor_occupancy, max_donor_occupancy = self.extract_donor_lone_pair_occupancy(self.min_donor_element, self.max_donor_element)
        return min_donor_occupancy, max_donor_occupancy

    def calculate_dipole_moment(self):
        dipole_vector = self.extract_dipole_moment()[1]
        dipole = np.sqrt(dipole_vector.dot(dipole_vector))  # unit: Debye dipole moment
        return dipole

    def calculate_dispersion_energy(self):
        dispersion_energy = self.extract_dispersion_energies()[-1]  # unit: eV
        return dispersion_energy

    def calculate_natural_charges(self):
        natural_charges = self.extract_natural_charges()
        metal_center_charge = natural_charges[self.metal_center_idx]
        max_donor_charge = natural_charges[self.max_donor_idx]
        min_donor_charge = natural_charges[self.min_donor_idx]
        return metal_center_charge, min_donor_charge, max_donor_charge

    def calculate_mulliken_charges(self):
        mulliken_charges = self.extract_mulliken_charges()
        metal_center_charge = mulliken_charges[self.metal_center_idx]
        max_donor_charge = mulliken_charges[self.max_donor_idx]
        min_donor_charge = mulliken_charges[self.min_donor_idx]
        return metal_center_charge, min_donor_charge, max_donor_charge

    def calculate_electronic_descriptors(self):
        homo_energy, lumo_energy, homo_lumo_gap = self.extract_homo_lumo_gap()
        hardness = 0.5*(lumo_energy - homo_energy)
        softness = 1./hardness
        electronegativity = (-0.5*(lumo_energy + homo_energy))
        electrophilicity = (electronegativity**2)/(2*hardness)
        return homo_energy, lumo_energy, homo_lumo_gap, hardness, softness, electronegativity, electrophilicity


if __name__ == '__main__':
    from molecular_graph import molecular_graph

    complexes_to_calc_descriptors = glob.glob(os.path.join(os.getcwd(), 'Workflow', '*.log'))
    for complex in complexes_to_calc_descriptors:
        # dft = DFTExtractor(complex, 39, 5, 28)
        elements, coordinates = read_cclib(complex)
        coordinates = coordinates[-1]
        # coordinates = np.array(coordinates).reshape(-1, len(coordinates), 3)
        # coordinates = coordinates[0]
        # # coordinates = np.array([np.array(coord) for coord in coordinates])
        # # print(coordinates)
        # elements = convert_elements(elements, output='symbols')
        elements = np.array(elements)
        geom_type = 'BD'
        ligand_atoms, bidentate = molecular_graph(elements=elements, coords=coordinates, geom=geom_type)
        dft = DFTExtractor(complex, bidentate[0] + 1, bidentate[1] + 1, bidentate[2] + 1, metal_adduct='nbd')
        print(dft.check_nbd_back_carbon())
        print(dft.get_hydrogens_bonded_to_carbon_back_nbd())
        # print(dft.extract_thermodynamic_descriptors())
        # print(dft.calculate_min_donor_metal_orbital_occupation(), dft.calculate_min_donor_metal_anti_orbital_occupation())
        # print(dft.calculate_max_donor_metal_orbital_occupation(), dft.calculate_max_donor_metal_anti_orbital_occupation())
        # print(dft.calculate_min_donor_other_orbital_occupation())
        # print(dft.calculate_max_donor_other_orbital_occupation())
        # print(dft.calculate_min_donor_other_anti_orbital_occupation())
        # print(dft.calculate_max_donor_other_anti_orbital_occupation())
        # print(dft.calculate_dipole_moment())
        # print(dft.calculate_donor_lone_pair_occupancy())
        # print(dft.calculate_dispersion_energy())
        # print(dft.calculate_natural_charges())
        # print(dft.calculate_mulliken_charges())
        # print(dft.calculate_electronic_descriptors())



