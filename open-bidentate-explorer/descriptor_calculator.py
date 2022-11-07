# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #
# ToDo: implement mark's workflow with Morfeus, but use ChemSpaX first and then the resulting MetalLigandComplex object

# from multiprocessing import freeze_support
# import datetime
# from rdkit import Chem
# from mordred import Calculator, descriptors
# import pandas as pd
#
#
# from model_structure_generator import MetalLigandComplex
#
#
# class DescriptorCalculator:
#     def __init__(self, ligand_dataframe):
#         self.calc = Calculator(descriptors, ignore_3D=True)
#         self.ligand_dataframe = ligand_dataframe
#
#     def calculate_mordred_descriptors(self, list_of_mols, list_of_identifier_columns):
#
#         # initialize 2D descriptor calculator
#         descriptors_df = self.calc.pandas(list_of_mols, nproc=1)
#
#         # enrich dataset
#         descriptors_df[list_of_identifier_columns] = self.ligand_dataframe[list_of_identifier_columns]
#         # add updatetime to data
#         descriptors_df['UpdateTime'] = datetime.datetime.now()
#
#         # write dataset with long errors (for human readability)
#         descriptors_df.to_csv(f'{}_long_errors.csv', encoding='utf-8')
#         # replace long errors with np.nan (for easier computational processing)
#         descriptors_df = descriptors_df.fill_missing()
#         descriptors_df.to_csv(f'{}.csv', encoding='utf-8')
#
#
#
#
#
#
#
# if __name__ == "__main__":
#     # needed for multiprocessing in Mordred to work correctly
#     freeze_support()
#
#     # initialize 2D descriptor calculator
#     calc = Calculator(descriptors, ignore_3D=True)
#
#     ligand_df = pd.read_excel("AH Catalyst Set and Substrate for MF.xlsx", "Sheet1")
#     smiles = ligand_df['smiles'].str.strip()
#
#     # calculate descriptors for each smiles and put the final result in a modified df
#     ligand_mols = [Chem.MolFromSmiles(smi) for smi in smiles]
#     ligand_descriptors_df = calc.pandas(ligand_mols)
#
#     # add name, aka, cas and smiles from original dataset
#     ligand_descriptors_df['name'] = ligand_df['Name']
#     ligand_descriptors_df['aka'] = ligand_df['Aka']
#     ligand_descriptors_df['cas'] = ligand_df['Cas']
#     ligand_descriptors_df['smiles'] = smiles
#
#     ligand_descriptors_df.to_csv('ligands_descriptors_long_errors.csv', encoding='utf-8')
#
#     # replace error messages with np.nan
#     ligand_descriptors_df = ligand_descriptors_df.fill_missing()
#     ligand_descriptors_df.to_csv('ligands_descriptors.csv', encoding='utf-8')
