# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #
# use the extracted xyz (via obelix) and the DFT log file of the free ligand to calculate descriptors
# generally only the electronic descriptors are extracted from the DFT log file since the steric/geometric descriptors
# are already calculated for the ligand (only the ligand) in the complex (via obelix)
import ast
import os

import numpy as np
import pandas as pd
import xtb.utils  # comment this when working on windows
import morfeus
from morfeus.io import write_xyz, read_cclib

from obelix.dft_extraction import DFTExtractor
from obelix.descriptor_calculator import Descriptors
from obelix.tools.utilities import dataframe_from_dictionary


class FreeLigand:
    def __init__(self, path_to_workflow, xyz_filename, dft_log_file):
        # the code expects the extracted xyz file of the free ligand and the DFT log file in the same directory
        self.path_to_workflow = path_to_workflow
        self.xyz_filename = xyz_filename

        # it is possible that the bidentate indices are different in the free ligand and complex, hence we need to
        # read the bidentate indices from the extraced ligand's xyz file (stored on the second line by obelix)
        self.complex_xyz_bidentate_1_idx = None
        self.complex_xyz_bidentate_2_idx = None
        self.free_ligand_xyz_bidentate_1_idx = None
        self.free_ligand_xyz_bidentate_2_idx = None

        if self.xyz_filename is not None:
            self.xyz_filename = os.path.join(path_to_workflow, xyz_filename)
            # read the string on the second line of the xyz file and try to convert to dict
            with open(self.xyz_filename, 'r') as f:
                line = f.readline()
                line = f.readline()
                try:
                    line = ast.literal_eval(line)
                    # read the 0-indexed bidentate indices from the dict and store in class
                    # we add + 1 to it since Morfeus and DFTExtractor expect 1-indexed indices
                    self.complex_xyz_bidentate_1_idx = line['complex_bidentate_1_index'] + 1
                    self.complex_xyz_bidentate_2_idx = line['complex_bidentate_2_index'] + 1
                    self.free_ligand_xyz_bidentate_1_idx = line['free_ligand_bidentate_1_index'] + 1
                    self.free_ligand_xyz_bidentate_2_idx = line['free_ligand_bidentate_2_index'] + 1
                except ValueError:
                    raise ValueError('The second line of the free ligand xyz file is not a dict')
        self.dft_log_file = os.path.join(path_to_workflow, dft_log_file)

        self.dft = None # the DFTExtractor class for this free ligand
        self.descriptor_calculator = Descriptors(central_atom=None,
                                                 path_to_workflow=self.path_to_workflow,
                                                 output_type='gaussian')
        # self.descriptor_df = None # the dataframe containing the descriptors for a set of free ligands

        # these will be used by morfeus, so we need to convert to 1-indexed later
        self.metal_center_idx = None  # this will remain None since we are not using a metal center
        self._min_donor_idx = None
        self._max_donor_idx = None

    def initialize_dft_extractor(self, log_file, metal_center_idx, min_donor_idx, max_donor_idx, metal_adduct):
        # initialize the DFTExtractor class using the given log file and the bidentate indices
        # having this function allows us to give the bidentate indices as input or determine them
        dft = DFTExtractor(log_file, metal_center_idx, min_donor_idx, max_donor_idx, metal_adduct)
        self.dft = dft
        return dft

    def assign_min_max_donor_xtb(self, elements, coordinates, calculation_method='gfn2_xtb', solvent=None):
        xtb = morfeus.XTB(elements, coordinates, solvent)
        atomic_charges = xtb.get_charges()  # determine min/max donor based on xTB charge of the donor atoms
        charge_bidentate_1 = atomic_charges[self.free_ligand_xyz_bidentate_1_idx]
        charge_bidentate_2 = atomic_charges[self.free_ligand_xyz_bidentate_2_idx]

        if charge_bidentate_1 > charge_bidentate_2:
            # this means that bidentate 1 is max donor (charge is less negative, so stronger donor)
            self._max_donor_idx = self.free_ligand_xyz_bidentate_1_idx
            self._min_donor_idx = self.free_ligand_xyz_bidentate_2_idx
        else:
            # this means that bidentate 2 is max donor
            self._max_donor_idx = self.free_ligand_xyz_bidentate_2_idx
            self._min_donor_idx = self.free_ligand_xyz_bidentate_1_idx

    def assign_min_max_donor_dft(self):
        if self.dft is None:
            raise ValueError('Initialize the DFTExtractor class first')
        # determine min/max donor based on DFT charge of the donor atoms
        natural_charges = self.dft.extract_natural_charges()
        charge_bidentate_1 = natural_charges[self.free_ligand_xyz_bidentate_1_idx]
        charge_bidentate_2 = natural_charges[self.free_ligand_xyz_bidentate_2_idx]

        if charge_bidentate_1 > charge_bidentate_2:
            # this means that bidentate 1 is max donor (charge is less negative, so stronger donor)
            self._max_donor_idx = self.free_ligand_xyz_bidentate_1_idx
            self._min_donor_idx = self.free_ligand_xyz_bidentate_2_idx
        else:
            # this means that bidentate 2 is max donor
            self._max_donor_idx = self.free_ligand_xyz_bidentate_2_idx
            self._min_donor_idx = self.free_ligand_xyz_bidentate_1_idx

        # reinitialize the DFTExtractor class with the new min/max donor indices
        self.initialize_dft_extractor(self.dft_log_file, self.metal_center_idx, self._min_donor_idx, self._max_donor_idx, metal_adduct='pristine')

    # it's possible that the min/max donor can be set externally (e.g. from an excel file with min/max from the complex)
    @property
    def min_donor_idx(self):
        return self._min_donor_idx

    @min_donor_idx.setter
    def min_donor_idx(self, min_donor_idx):
        # ToDo: check if this is a valid index before setting it
        self._min_donor_idx = min_donor_idx

    @property
    def max_donor_idx(self):
        return self._max_donor_idx

    @max_donor_idx.setter
    def max_donor_idx(self, max_donor_idx):
        # ToDo: check if this is a valid index before setting it
        self._max_donor_idx = max_donor_idx

    # use the Descriptor class and _calculate_dft_descriptors_from_log to calculate the DFT descriptors
    def calculate_dft_descriptors(self, solvent=None, extract_xyz_from_log=False, printout=False):
        # for the free ligand we need to have the extracted xyz file to get the bidentate indices and the matching
        # log file to get the DFT data, hence we cannot iterate over log files like we did in obelix.Descriptors.calculate_dft_descriptors
        properties = {}
        base_with_extension = os.path.basename(self.dft_log_file)
        split_base = os.path.splitext(base_with_extension)
        filename = split_base[0]
        print('\nCalculating descriptors for:', filename)
        properties['filename_tud'] = filename
        properties['index_donor_min'] = self._min_donor_idx
        properties['index_donor_max'] = self._max_donor_idx

        elements, coordinates = None, None
        # error catching for reading elements and coordinates from log file
        try:
            elements, coordinates = read_cclib(self.dft_log_file)
            if not len(coordinates[-1]) == 3:  # if this is true, there is only 1 coordinates array
                coordinates = coordinates[
                    -1]  # else morfeus descriptors are calculated for last geometry in log file
            elements = np.array(elements)
        except:
            print('Error reading elements and coordinates from log file for: ', filename)
            print('Make sure to check the geometry')

        # calculate DFT descriptors from Gaussian log file
        # get indices of bidentate ligands and metal for descriptor calculation class
        dft_properties = {}
        try:
            if self.dft is None:
                print('Initializing DFTExtractor class')
                dft = self.initialize_dft_extractor(self.dft_log_file, self.metal_center_idx, self._min_donor_idx,
                                                    self._max_donor_idx, metal_adduct='pristine')
            else:
                dft = self.dft
            dft_properties = self.descriptor_calculator._calculate_dft_descriptors_from_log(dft, dft_properties)
        except Exception as e:
            print(e)
            print(f'DFT descriptor calculation failed for {filename}')

        if len(dft_properties) > 0:
            properties.update(dft_properties)

        # write xyz for log file
        if extract_xyz_from_log:
            dft_xyz_filename = self.dft_log_file[:-4] + '_DFT.xyz'
            if elements is not None and coordinates is not None:
                write_xyz(os.path.join(self.path_to_workflow, dft_xyz_filename), elements, coordinates)
            elif elements is None or coordinates is None:
                print('Error writing xyz for: ', filename)
                print('Make sure to check the geometry')

        # now all properties need to have free_ligand_ prefix to distinguish them from the complex descriptors
        for key in list(properties.keys()):
            properties[f'free_ligand_{key}'] = properties.pop(key)

        # merging all properties into one dictionary can now be done outside of this class
        # unlike in obelix.Descriptors.calculate_dft_descriptors
        return properties


if __name__ == "__main__":
    # example usage
    ligand_number_list = [f'L{i}' for i in range(1, 3)]
    dictionary_for_properties = {}

    for ligand_number in ligand_number_list:
        # files are located in ../tests/files
        path_to_test_files = os.path.join(os.getcwd(), '..', 'tests', 'files')
        free_ligand = FreeLigand(path_to_test_files, f'free_ligand_{ligand_number}.xyz', f'free_ligand_{ligand_number}_SP.log')
        free_ligand.initialize_dft_extractor(free_ligand.dft_log_file, None, free_ligand._min_donor_idx, free_ligand._max_donor_idx, metal_adduct='pristine')
        # print(free_ligand.min_donor_idx)
        free_ligand.assign_min_max_donor_dft()
        # print(free_ligand.min_donor_idx)
        properties = free_ligand.calculate_dft_descriptors()
        dictionary_for_properties[ligand_number] = properties
    new_descriptor_df = dataframe_from_dictionary(dictionary_for_properties)
    # reset the index and name that column to 'Ligand#' for consistency with the complex descriptor df
    new_descriptor_df = new_descriptor_df.reset_index().rename(columns={'index': 'Ligand#'})
    new_descriptor_df.to_csv('test_free_ligand_descriptors.csv', index=False)
