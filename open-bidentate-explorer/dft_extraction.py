# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #
import os
import glob
from morfeus.io import read_cclib
from cclib.parser import ccopen


class DFTExtractor(object):
    def __init__(self, log_file, metal_center_idx, min_donor_idx, max_donor_idx):
        self.log_file = log_file
        # id's come from morfeus, so they start at 1, subtract 1 to get the correct index for cclib
        self.metal_center_idx = metal_center_idx - 1
        self.min_donor_idx = min_donor_idx - 1
        self.max_donor_idx = max_donor_idx - 1

        # use cclib parser
        self.parser = ccopen(log_file)
        self.data = self.parser.parse()
        self.atom_charges_dict = self.data.atomcharges

        # use morfeus parser
        self.elements, self.coordinates = read_cclib(log_file)

        # read the log file for extraction
        with open(log_file) as f:
            self.log_file_text = f.readlines()

    def extract_homo_lumo_gap(self):
        homo = self.data.homos[0]
        energies = self.data.moenergies[0]
        lumo_energy = energies[homo+1]
        homo_energy = energies[homo]
        homo_lumo_gap = lumo_energy - homo_energy
        return homo_energy, lumo_energy, homo_lumo_gap

    def extract_natural_charges(self):
        natural_charges = self.atom_charges_dict['natural']
        return natural_charges

    def extract_mulliken_charges(self):
        atomic_charges = self.atom_charges_dict['mulliken']
        return atomic_charges

    def extract_energy(self):
        pass

    def extract_thermodynamic_descriptors(self):
        # extract energy, enthalpy, entropy, free energy
        free_energy = self.data.freeenergy
        enthalpy = self.data.enthalpy
        zero_point_correction = self.data.zpve
        entropy = self.data.entropy
        pass

    def calculate_natural_charges(self):
        natural_charges = self.extract_natural_charges()
        metal_center_charge = natural_charges[self.metal_center_idx]
        max_donor_charge = natural_charges[self.max_donor_idx]
        min_donor_charge = natural_charges[self.min_donor_idx]
        return metal_center_charge, max_donor_charge, min_donor_charge

    def calculate_mulliken_charges(self):
        mulliken_charges = self.extract_mulliken_charges()
        metal_center_charge = mulliken_charges[self.metal_center_idx]
        max_donor_charge = mulliken_charges[self.max_donor_idx]
        min_donor_charge = mulliken_charges[self.min_donor_idx]
        return metal_center_charge, max_donor_charge, min_donor_charge

    def calculate_electronic_descriptors(self):
        homo_energy, lumo_energy, homo_lumo_gap = self.extract_homo_lumo_gap()
        hardness = 0.5*(lumo_energy - homo_energy)
        softness = 1./hardness
        electronegativity = (-0.5*(lumo_energy + homo_energy))
        electrophilicity = (electronegativity**2)/(2*hardness)
        return homo_energy, lumo_energy, homo_lumo_gap, hardness, softness, electronegativity, electrophilicity


if __name__ == '__main__':
    complexes_to_calc_descriptors = glob.glob(os.path.join(os.getcwd(), 'Workflow', '*.log'))
    for complex in complexes_to_calc_descriptors:
        dft = DFTExtractor(complex, 1, 1, 1)
        # print(dft.extract_homo_lumo_gap())
        print(dft.extract_natural_charges())
