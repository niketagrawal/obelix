# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #
import os
import glob
from morfeus.io import read_cclib
import cclib


class DFTExtractor(object):
    def __init__(self, log_file):
        self.log_file = log_file

    def extract_energy(self):
        pass

    def extract_thermodynamic_descriptors(self):
        # extract energy, enthalpy, entropy, free energy
        pass

    def extract_electronic_descriptors(self):
        pass


if __name__ == '__main__':
    complexes_to_calc_descriptors = glob.glob(os.path.join(os.getcwd(), 'Workflow', '*.log'))
    for complex in complexes_to_calc_descriptors:
        elements, coordinates = read_cclib(complex)
