# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #
import mace
from rdkit import Chem
import pandas as pd
import os


class MetalLigandComplex:
    def __init__(self, ligand_dataframe, substrate, central_atom):
        self.ligand_df = ligand_dataframe
        self.substrate = substrate
        self.central_atom = central_atom

    @staticmethod
    def make_dir(folder_name):
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

    @staticmethod
    def transform_smiles_to_mace_input(list_of_ligands):
        # smarts is used in substructure search
        # donor_atoms_smarts = ['[PX3]', '[NX3]']  # search for P or N with 3 neighbours
        donor_atoms_smarts = ['[PX3]']  # search for P with 3 neighbours

        list_of_mace_ligands = []
        for i in list_of_ligands:
            mol = Chem.MolFromSmiles(i)
            for smarts in donor_atoms_smarts:
                # get the donor atom index from the tuple returned by GetSubstructureMatches
                donor_atom_idxs = [idx[0] for idx in mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))]
                for idx in donor_atom_idxs:
                    # set atom map number in SMILES, example: PCCP --> [P:1]CC[P:1]
                    mol.GetAtomWithIdx(idx).SetAtomMapNum(1)
                list_of_mace_ligands.append(Chem.MolToSmiles(mol))

        mace_ligands = list(set(list_of_mace_ligands))
        return mace_ligands

    def create_complex_mace(self, identifier_column, ligand_smiles_column, auxillary_ligands, CA,
                            geom, output_file_identifier):
        bidentate_ligands = self.ligand_df[ligand_smiles_column].dropna()
        # bidentate_ligands = self.transform_smiles_to_mace_input(bidentate_ligands)
        list_of_identifier = self.ligand_df[identifier_column]

        for index, ligand in enumerate(bidentate_ligands):
            ligands = [ligand]
            ligands.extend(auxillary_ligands)

            X = mace.ComplexFromLigands(ligands, CA, geom)
            Xs = X.GetStereomers(regime='all', dropEnantiomers=True)
            print(f'OH Stereomers for ligand {list_of_identifier[index]}:', len(Xs))

            for i, X in enumerate(Xs):
                X.AddConformers(numConfs=10)
                self.make_dir(f'ligands/{list_of_identifier[index]}')
                X.ToXYZ(f'ligands/{list_of_identifier[index]}/{output_file_identifier}_{self.central_atom}_{i}.xyz',
                        confId='min')

    def create_octahedral_structures(self, identifier_column, ligand_smiles_column):
        list_of_auxillary_ligands = ['[H-:1]', '[H-:1]', '[H-:1]', self.substrate]
        print("Octahedral structure with substrate")
        output_file_identifier = 'OH'
        self.create_complex_mace(identifier_column=identifier_column, ligand_smiles_column=ligand_smiles_column,
                                 auxillary_ligands=list_of_auxillary_ligands, CA=self.central_atom, geom='OH',
                                 output_file_identifier=output_file_identifier)

    def create_octahedral_structures_without_substrate(self, identifier_column, ligand_smiles_column):
        # needed if you want to calculate energy of reaction
        print("Octahedral structure without substrate")
        list_of_auxillary_ligands = ['[H-:1]', '[H-:1]', '[H-:1]']
        output_file_identifier = 'OH_no_substrate'
        self.create_complex_mace(identifier_column=identifier_column, ligand_smiles_column=ligand_smiles_column,
                                 auxillary_ligands=list_of_auxillary_ligands, CA=self.central_atom, geom='OH',
                                 output_file_identifier=output_file_identifier)

    def create_square_planar_structures_without_substrate(self, identifier_column, ligand_smiles_column):
        print("SP without substrate")
        list_of_auxillary_ligands = []
        output_file_identifier = 'SP_no_substrate'
        self.create_complex_mace(identifier_column=identifier_column, ligand_smiles_column=ligand_smiles_column,
                                 auxillary_ligands=list_of_auxillary_ligands, CA=self.central_atom, geom='SP',
                                 output_file_identifier=output_file_identifier)

    def create_square_planar_structures(self, identifier_column, ligand_smiles_column):
        print("SP with substrate")
        list_of_auxillary_ligands = [self.substrate]
        output_file_identifier = 'SP'
        self.create_complex_mace(identifier_column=identifier_column, ligand_smiles_column=ligand_smiles_column,
                                 auxillary_ligands=list_of_auxillary_ligands, CA=self.central_atom, geom='SP',
                                 output_file_identifier=output_file_identifier)


if __name__ == "__main__":
    dataset_name = "catalyst_set.xlsx"
    tab_name = "selected_ligands"
    column_name = "smiles_for_mace"

    dataset_df = pd.read_excel(dataset_name, tab_name)
    imine_substrate = 'CC#[N:1]'  # acetonnitril
    central_atom = '[Ir+3]'
    mace_complex = MetalLigandComplex(ligand_dataframe=dataset_df, substrate=imine_substrate,
                                      central_atom=central_atom)

    # The four functions above (SP and OH with and without structure can be called through the 4 lines below)
    mace_complex.create_octahedral_structures(identifier_column='index', ligand_smiles_column=column_name)
    # mace_complex.create_octahedral_structures_without_substrate(identifier_column='index', ligand_smiles_column=columnName)
    # mace_complex.create_square_planar_structures(identifier_column='index', ligand_smiles_column=columnName)
    mace_complex.create_square_planar_structures_without_substrate(identifier_column='index',
                                                                   ligand_smiles_column=column_name)

