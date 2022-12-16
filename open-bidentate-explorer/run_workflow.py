from chemspax.main import main
from chemspax.utilities import *
import os
import glob
import pandas as pd
from tools.utilities import *
import shutil
from tqdm import tqdm
import mace
import subprocess
import morfeus as mf
from morfeus.conformer import ConformerEnsemble
from morfeus import BiteAngle, ConeAngle, BuriedVolume, Dispersion, SASA, read_xyz
from xtb.utils import get_method
import pandas as pd
from openbabel import openbabel
from molecular_graph import molecular_graph
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.info')


class MACE:   
    def __init__(self, bidentate, CA, name_of_xyz):
        # Declare instances of the class (a.k.a self variable)
        self.CA = CA
        self.name_of_xyz = name_of_xyz
        self.bidentate = bidentate

    @staticmethod
    def find_bidentate(mol, metal_atom_atomic_number):
        """Find indices of bidentate ligand if functions implemented in MACE class are not working or if
        structure was made manually and there is no mapped SMILES string

        Args:
            mol: RDKit mol object
            metal_atom_atomic_number:

        Returns: list of bidentate atom indices

        """
        bidentate_idxs = []
        try:
            # find bidentate indices for first conformer, they should be the same for all conformers
            # in case of SP we only need to find metal atom and neighbours since only one cycle with metal center is possible
            metal_atom = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == metal_atom_atomic_number][0]
            metal_bonds = [bond for bond in metal_atom.GetBonds()]
            for bond in metal_bonds:
                # atoms that are connected to the metal and in the bidentate cycle are the bidentate atoms
                idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                bidentate_idxs.append(idx1) if idx1 != metal_atom.GetIdx() else bidentate_idxs.append(idx2)

        except Exception as e:
            # ToDo: add logging
            # go on to next conformer and try again
            return None

        return bidentate_idxs[0], metal_atom.GetIdX(), bidentate_idxs[1]

    def generate_complex_SP_xyz(self):
        geom = 'SP'

        core = mace.ComplexFromLigands([self.bidentate], self.CA, geom)
        Xs = core.GetStereomers(regime='all', dropEnantiomers=True)

        bidentate_idxs = []  # the final bidentate cycle indices
        for i, X in enumerate(Xs):
            X.AddConformers(numConfs=10)   
            X.ToXYZ(self.CA + '_' + '{}{}.xyz'.format(self.name_of_xyz, i), confId='min')
            if len(bidentate_idxs) == 0:
                try:
                    # find bidentate indices for first conformer, they should be the same for all conformers
                    # in case of SP we only need to find metal atom and neighbours since only one cycle with metal center is possible
                    atomic_number_metal_center = Chem.MolFromSmiles(self.CA).GetAtomWithIdx(0).GetAtomicNum()
                    metal_atom = [atom for atom in X.mol3D.GetAtoms() if atom.GetAtomicNum() == atomic_number_metal_center][0]
                    metal_bonds = [bond for bond in metal_atom.GetBonds()]
                    for bond in metal_bonds:
                        # atoms that are connected to the metal and in the bidentate cycle are the bidentate atoms
                        idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                        bidentate_idxs.append(idx1) if idx1 != metal_atom.GetIdx() else bidentate_idxs.append(idx2)

                except Exception as e:
                    # ToDo: add logging
                    # go on to next conformer and try again
                    continue

        return bidentate_idxs[0], metal_atom.GetIdX(), bidentate_idxs[1]

    def generate_complex_OH_xyz(self, auxiliary_ligands = [], substrate = []):
        geom = 'OH'

        # Check if auxiliary ligands and substrate are present
        if auxiliary_ligands == [] and substrate == []:
            auxiliary_ligands = ['[H-:1]']*4
        if auxiliary_ligands == [] and substrate != []:
            auxiliary_ligands = ['[H-:1]']*3
        
        auxiliary = auxiliary_ligands
        
        if substrate == []:
            auxiliary.extend([self.bidentate])
        else:
            auxiliary.extend([self.bidentate, substrate[0]])

        core = mace.ComplexFromLigands(auxiliary, self.CA, geom)
        Xs = core.GetStereomers(regime='all', dropEnantiomers=True)

        bidentate_cycle_idxs = None  # all cycles that could possibly be the metal-bidentate cycle
        bidentate_idxs = []  # the final bidentate cycle indices
        for i, X in enumerate(Xs):
            X.AddConformers(numConfs=10)   
            X.ToXYZ(self.CA + '_' + '{}{}.xyz'.format(self.name_of_xyz, i), confId='min')

            # find bidentate indices on the conformer if none were yet
            if len(bidentate_idxs) == 0:
                try:
                    # get all mapped (donor) atom indices
                    donor_atoms = [atom.GetIdx() for atom in X.mol3D.GetAtoms() if atom.GetAtomMapNum()]
                    # smallest set of simple rings containing at least two mapped atoms
                    # Chem.GetSSSR() is updated in rdkit 2020.09.1, versions before had to use Chem.GetSymmSSSR()
                    atoms_of_smallest_rings = [list(idx) for idx in Chem.GetSSSR(X.mol3D)]
                    for list_idxs in atoms_of_smallest_rings:
                        # if at least two mapped atoms are in the ring it's a multidentate ligand
                        if check_if_at_least_two_mapped_atoms_in_ring(donor_atoms, list_idxs):
                            bidentate_cycle_idxs = list_idxs
                            break

                    # find metal atom and neighbours
                    atomic_number_metal_center = Chem.MolFromSmiles(self.CA).GetAtomWithIdx(0).GetAtomicNum()
                    metal_atom = [atom for atom in X.mol3D.GetAtoms() if atom.GetAtomicNum() == atomic_number_metal_center][0]
                    metal_bonds = [bond for bond in metal_atom.GetBonds()]
                    for bond in metal_bonds:
                        # atoms that are connected to the metal and in the bidentate cycle are the bidentate atoms
                        idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                        if idx1 == metal_atom.GetIdx() and idx2 in bidentate_cycle_idxs:
                            bidentate_idxs.append(idx2)
                        elif idx1 in bidentate_cycle_idxs and idx2 == metal_atom.GetIdx():
                            bidentate_idxs.append(idx1)

                except Exception as e:
                    # ToDo: add logging
                    # go on to next conformer and try again
                    continue

        return bidentate_idxs[0], metal_atom.GetIdX(), bidentate_idxs[1]
            

class Workflow:
    def __init__(self, mace_input = [], chemspax_input = [], crest_input = [], path_to_workflow = []):
        
        self.mace_input = mace_input
        self.chemspax_input = chemspax_input
        self.crest_input = crest_input
        self.path_to_workflow = path_to_workflow

        self.bidentate_1_index = None
        self.bidentate_2_index = None
        self.metal_index = None
        
        # Unpack inputs of MACE, Chemspax, CREST
        print('Workflow is initializing. Converting your dict. input to variables.')
        print('')
        
        if mace_input != []:  
          self.mace_ligands, self.auxiliary_ligands, self.geom, self.central_atom, self.names_of_xyz, self.substrate = self.initialize_mace()
        if chemspax_input != []:        
          self.substituent_list, self.path_to_database, self.path_to_substituents = self.initialize_chemspax()
        if crest_input != []:
          self.method, self.charge_of_complex, self.multiplicity, self.solvent, self.conf_search = self.initialize_crest()
              
    def initialize_mace(self):
        print('Reading MACE inputs')
        
        bidentate = list(pd.read_excel(self.mace_input['bidentate_ligands'])['smiles'])
        auxiliary_ligands = self.mace_input['auxiliary_ligands']
        geom = self.mace_input['geom']
        central_atom = self.mace_input['central_atom']
        names_of_xyz_key = list(pd.read_excel(self.mace_input['bidentate_ligands'])['Name'])
        substrate = self.mace_input['substrate']
        return bidentate, auxiliary_ligands, geom, central_atom, names_of_xyz_key, substrate
      
    def initialize_chemspax(self):

        print('Reading ChemSpaX inputs')

        substituent_list = self.chemspax_input['substituent_list']
        path_to_database = self.chemspax_input['path_to_database']
        path_to_substituents = self.chemspax_input['path_to_substituents']
        
        return substituent_list, path_to_database, path_to_substituents
    
    def initialize_crest(self): 
        print('Reading CREST inputs.')
        
        method = self.crest_input['method']
        multiplicity = self.crest_input['multiplicity']
        charge_of_complex = self.crest_input['charge_complex']
        solvent = self.crest_input['solvent']
        conf_search = self.crest_input['conformer_search']
        return method, charge_of_complex, multiplicity, solvent, conf_search
        
    def prepare_folder_structure(self):
        
        print('Preparing the folder structure')
    
        os.mkdir(self.path_to_workflow)
        os.chdir(self.path_to_workflow)
        os.mkdir('MACE')
        os.mkdir('ChemSpaX')
        os.mkdir('CREST')
        os.mkdir('Descriptors')      
        os.chdir(self.path_to_workflow)

    def run_mace(self):
        
        os.chdir(os.path.join(self.path_to_workflow, 'MACE'))

        for i in tqdm(range(len(self.mace_ligands))):
            try:
              if self.geom == 'OH':
                  identifier_OH = self.names_of_xyz[i] + '_OH'
                  complex = MACE(self.mace_ligands[i], central_atom,  identifier_OH)
                  self.bidentate_1_index, self.metal_index, self.bidentate_2_index = \
                      complex.generate_complex_OH_xyz(auxiliary_ligands, substrate=self.substrate)
              
              if self.geom == 'SP':
                  identifier_SP = self.names_of_xyz[i] + '_SP'
                  complex = MACE(self.mace_ligands[i], central_atom,  identifier_SP)
                  self.bidentate_1_index, self.metal_index, self.bidentate_2_index = complex.generate_complex_SP_xyz()
            except Exception:
              print('Wrong SMILES:', i)               
    
    def run_chemspax(self, names, skeleton_list_):  
        
        """
        
        Run chemspax from the pip installable package chemspax.
        Use prepare_data(...) to prepare the substituents.
        Then run main(...) to functionalize the skeletons prepared by MACE.
        
        """
        
        print(self.substituent_list)
        chemspax_working_directory = os.path.join(self.path_to_workflow, 'ChemSpaX')
        
        # Export MACE skeletons to ChemSpaX
        
        src_dir_skeletons =  os.path.join(self.path_to_workflow, 'MACE')
        
        if self.substituent_list == []:
            dest_dir_skeletons = os.path.join(self.path_to_workflow, 'ChemSpaX', 'chemspax_output')
            shutil.copytree(src_dir_skeletons, dest_dir_skeletons)  

        else:              
            dest_dir_skeletons =  os.path.join(self.path_to_workflow, 'ChemSpaX', 'skeletons')        
            shutil.copytree(src_dir_skeletons, dest_dir_skeletons)
            skeleton_list = glob.glob(os.path.join(self.path_to_workflow , 'ChemSpaX', 'skeletons','*.xyz'))
        
            for xyz_file in skeleton_list:
                change_second_line_xyz(xyz_file, new_content='')

            # set skeleton path to ../ChemSpaX/skeletons        
            path_to_skeletons = os.path.join(os.getcwd(), dest_dir_skeletons)
            
            # copy substituents from user source to ChemSpaX folder
            src_dir_subs = self.path_to_substituents
            dest_dir_subs = os.path.join(self.path_to_workflow, 'ChemSpaX', 'substituents_xyz')        
            shutil.copytree(src_dir_subs, dest_dir_subs)

            # Prepare data (from data_preparation.py in the chemspax package)
            print('Data preparation has been performed.')
            path_to_substituents = os.path.join(os.getcwd(), dest_dir_subs)
            #   prepare_data(path_to_substituents, path_to_skeletons, self.path_to_database)

            # Set path to output
            path_to_output = os.path.join(self.path_to_workflow, 'ChemSpaX','chemspax_output')
        
            # Run chemspax from main
        for index, sub_list in enumerate(self.substituent_list):
            #### Run chemspa
            main([os.path.join(path_to_skeletons, skeleton_list_[index])], sub_list, self.path_to_database, path_to_substituents, os.path.join(path_to_skeletons), chemspax_working_directory, path_to_output)
            
            os.rename(os.path.join(path_to_output, skeleton_list_[index][:-4] + '_func_' + str(len(sub_list)) + '.mol'), \
                os.path.join(path_to_output, list(names)[index] +  '.mol'))
            
            ## Convert mol to xyz to keep the bonding information
            obconversion = openbabel.OBConversion()
            obconversion.SetInFormat('mol')
            obconversion.SetOutFormat('xyz')
            mol = openbabel.OBMol()
            obconversion.ReadFile(mol, \
                os.path.join(path_to_output, list(names)[index] +  '.mol'))
            
            obconversion.WriteFile(mol, \
                os.path.join(path_to_output, list(names)[index] +  '.xyz'))
            
            for filename in glob.glob(os.path.join(self.path_to_workflow, 'ChemSpaX', 'chemspax_output', '*_func_*')):
                os.remove(filename) 
    
    def run_crest(self, path_to_complexes = [], path_to_output = [], conformer_search = 'off'):
        
        """
        
        This function is designed to either run conformer search with xtb preopt. or just xtb optimization.
        The subprocess.call is used to make a binary call from your ./bashrc file or any local path in your 
        $PATH system variable.
        
        Files are stored in the users preferred place when you run the workflow in a modular way, or by design
        the results are stored in your Workflow/CREST folder.

        """
    
        
        if path_to_complexes == []:
          path_to_complexes = os.path.join(self.path_to_workflow, 'ChemSpaX', 'chemspax_output')
          os.chdir(self.path_to_workflow)

        
        for path_to_complex in glob.glob(os.path.join(path_to_complexes, '*.xyz')):
            os.chdir(path_to_complexes)
            
            folder_name = path_to_complex[:-4]
            os.mkdir(folder_name)
            os.chdir(folder_name)

            print('About to run {}-xtb optimization:'.self.method)

            if conformer_search == 'on':
                subprocess.call('xtb "%s" --%s --chrg %s  --uhf %s --alpb %s --opt > xtb.out' \
                    % (path_to_complex, self.method, self.charge_of_complex, self.multiplicity, self.solvent), shell=True)
                
                print('About to RUN conformer search with CREST:')
                
                subprocess.call('crest xtbopt.xyz --%s --chrg %s  --uhf %s --alpb %s > crest.out' \
                    % (path_to_complex, self.method, self.charge_of_complex, self.multiplicity, self.solvent), shell=True)
       
            else:
                subprocess.call('xtb "%s" --%s --chrg %s  --uhf %s --alpb %s --opt > xtb.out' \
                    % (path_to_complex, self.method, self.charge_of_complex, self.multiplicity, self.solvent), shell=True)
            
            
            src_dir = path_to_complex[:-4]
            
            if path_to_output != []:
                dest_dir = path_to_output
            else:
                dest_dir = os.path.join(self.path_to_workflow, 'CREST')
            
            shutil.move(src_dir, dest_dir) 
            os.chdir(path_to_complexes)
        
        os.chdir(self.path_to_workflow)
            
    def calculate_descriptors(self, path_to_outputs = [], output = 'xtb'):
        
        """
        
        This function calculates morfeus descriptors.
        
        """
        
        if path_to_outputs == []:
            complexes_to_calc_descriptors = glob.glob(os.path.join(self.path_to_workflow, 'CREST', '*'))
        else:
            complexes_to_calc_descriptors = glob.glob(os.path.join(path_to_outputs, '*'))
        # Check output     
        ### Data structure 'complex': {}
        
        # Get descriptors from xtb-optimized structures
        if output == 'xtb':
            dictionary_for_properties = {}
            
            for complex in complexes_to_calc_descriptors:
                properties = {}
                print("Descriptor calculation failed for this complex:", os.path.basename(os.path.normpath(complex)))   

                elements, coordinates = read_xyz(os.path.join(complex, 'xtbopt.xyz'))
                print(len(elements))
                # if self.bidentate_1_index is None or self.bidentate_2_index is None or self.metal_index is None:
                #     bidentate = find_bidentate(os.path.join(complex, 'xtbopt.xyz'))
                # else:
                #     bidentate = [self.bidentate_1_index, self.metal_index, self.bidentate_2_index]
                
                
                ligand_atoms, bidentate = molecular_graph(elements=elements, coords=coordinates)
                #     # ToDo: use find_bidentate static method from MACE class if None, but need mol object returned
                
                #### Set up a as bidentate + metal array
                a = [bidentate[0]]
                a.extend(ligand_atoms[bidentate[1]])
                a = list(np.sort(np.array(a)))
                
                ### Get index of metal in the structure without the aux ligands and the substrate                
                for id, i in enumerate(a):
                    if i == bidentate[0]:
                        diff = id
                                 
                # Diff is the new index of the Metal 
                # diff + 1 from the indexing in morfeus ( + 1)
                                       
                properties["bite_angle"] = BiteAngle(coordinates, bidentate[0] + 1, bidentate[1] + 1, bidentate[2] + 1).angle
                elements_cone_angle = elements[a]                
                coordinates_cone_angle = np.array(coordinates)[a]
                
                properties["cone_angle"] = ConeAngle(elements_cone_angle, coordinates_cone_angle, diff + 1).cone_angle
                
                bv1 = BuriedVolume(elements, coordinates, bidentate[0] + 1, radius=3.5)
                bv2 = BuriedVolume(elements, coordinates, bidentate[1] + 1, radius=3.5)
                bv3 = BuriedVolume(elements, coordinates, bidentate[2] + 1, radius=3.5)
                
                
                properties["buried_volume_metal_center"] = bv1.fraction_buried_volume
                properties["buried_volume_P1"] = bv2.fraction_buried_volume
                properties["buried_volume_P2"] = bv3.fraction_buried_volume

                # print(BuriedVolume(ce.elements, conformer.coordinates, bidentate[1]).print_report())
                # BuriedVolume(elements, coordinates, bidentate[1], radius=4).print_report()

                properties["dispersion"] = Dispersion(elements, coordinates).atom_p_int[bidentate[0] + 1]
                properties["sasa"] = SASA(elements, coordinates).area
                
                # Electronic descriptors
                
                instance_electronic = mf.XTB(elements, coordinates, solvent='ch2cl2')
                
                properties["ip"] = instance_electronic.get_ip()
                properties["dipole"] = instance_electronic.get_dipole().dot(instance_electronic.get_dipole())
                properties["ea"] = instance_electronic.get_ea()     
                properties["electrofugality"] = instance_electronic.get_global_descriptor(variety = 'electrofugality')
                properties["nucleofugality"] = instance_electronic.get_global_descriptor(variety = 'nucleofugality')
                properties["nucleophilicity"] = instance_electronic.get_global_descriptor(variety = 'nucleophilicity')
                properties["electrophilicity"] = instance_electronic.get_global_descriptor(variety = 'electrophilicity')

                homo = instance_electronic.get_homo()
                lumo = instance_electronic.get_lumo()    
                properties["HOMO_LUMO_gap"] = lumo - homo
                
                for property in properties.keys():
                    dictionary_for_properties[os.path.basename(os.path.normpath(complex))] = properties

            # except Exception:
            #     print("Descriptor calculation failed for this complex:", os.path.basename(os.path.normpath(complex)))   
                         
        # Iterate through the CREST outputs of different descriptors
        else:
            dictionary_for_properties = {}    
            for complex in complexes_to_calc_descriptors:
                ce = ConformerEnsemble.from_crest(complex)
                ce.prune_rmsd()
                ce.sort()
                for conformer in ce:        
                    # sasa = ce.boltzmann_statistic("sasa")
                    # sasa_std = ce.boltzmann_statistic("sasa", statistic = "std")

                    # if self.bidentate_1_index is None or self.bidentate_2_index is None or self.metal_index is None:
                    #     bidentate = find_bidentate(conformer)
                    # else:
                    #     bidentate = [self.bidentate_1_index, self.metal_index, self.bidentate_2_index]
                    # #     # ToDo: use find_bidentate static method from MACE class if None, but need mol object returned

                                    #### Set up a as bidentate + metal array
                    a = [bidentate[0]]
                    a.extend(ligand_atoms[bidentate[1]])
                    a = list(np.sort(np.array(a)))
                    
                    ### Get index of metal in the structure without the aux ligands and the substrate                
                    for id, i in enumerate(a):
                        if i == bidentate[0]:
                            diff = id
                                    
                    # Diff is the new index of the Metal 
                    # diff + 1 from the indexing in morfeus ( + 1)
                    properties["bite_angle"] = BiteAngle(coordinates, bidentate[0] + 1, bidentate[1] + 1, bidentate[2] + 1).angle
                    elements_cone_angle = elements[a]                
                    coordinates_cone_angle = np.array(coordinates)[a]
                
                    # bidentate = find_bidentate(conformer)
           
                    conformer.properties["bite_angle"] = BiteAngle(conformer.coordinates, bidentate[1], bidentate[0], bidentate[2]).angle
                    conformer.properties["cone_angle"] = ConeAngle(elements_no_H, coordinates_no_H, bidentate[1]).cone_angle
                    conformer.properties["buried_volume"] = BuriedVolume(ce.elements, conformer.coordinates, bidentate[1]).fraction_buried_volume
                    
                    conformer.properties["cone_angle"] = ConeAngle(elements_cone_angle, coordinates_cone_angle, diff + 1).cone_angle
                
                    bv1 = BuriedVolume(elements, coordinates, bidentate[0] + 1, radius=3.5)
                    bv2 = BuriedVolume(elements, coordinates, bidentate[1] + 1, radius=3.5)
                    bv3 = BuriedVolume(elements, coordinates, bidentate[2] + 1, radius=3.5)
                    
                    
                    properties["buried_volume_metal_center"] = bv1.fraction_buried_volume
                    properties["buried_volume_P1"] = bv2.fraction_buried_volume
                    properties["buried_volume_P2"] = bv3.fraction_buried_volume


                    properties["dispersion"] = Dispersion(elements, coordinates).atom_p_int[bidentate[0] + 1]
                    properties["sasa"] = SASA(elements, coordinates).area
                    # Electronic descriptors
                    
                    instance_electronic = mf.XTB(ce.elements, conformer.coordinates, solvent='ch2cl2')
                    conformer.properties["ip"] = instance_electronic.get_ip()
                    conformer.properties["dipole"] = instance_electronic.get_dipole().dot(instance_electronic.get_dipole())
                    conformer.properties["ea"] = instance_electronic.get_ea()
                    conformer.properties["electrofugality"] = instance_electronic.get_global_descriptor(variety = 'electrofugality')
                    conformer.properties["nucleofugality"] = instance_electronic.get_global_descriptor(variety = 'nucleofugality')
                    conformer.properties["nucleophilicity"] = instance_electronic.get_global_descriptor(variety = 'nucleophilicity')
                    conformer.properties["electrophilicity"] = instance_electronic.get_global_descriptor(variety = 'electrophilicity')
                    
                    homo = instance_electronic.get_homo()
                    lumo = instance_electronic.get_lumo()    
                    conformer.properties["HOMO_LUMO_gap"] = homo - lumo

                for property in ce.get_properties().keys():
                    dictionary_for_properties[os.path.basename(os.path.normpath(complex))] = {property: ce.boltzmann_statistic(property)}
                
        dataframe = dataframe_from_dictionary(dictionary_for_properties)
        dataframe.to_excel('descriptors_OH_oct_quads.xlsx')     
        
    def run_workflow(self):
        
        if self.path_to_workflow != []:
           self.prepare_folder_structure()
        
        self.run_mace()
        self.run_chemspax()
        self.run_crest(self.conf_search)
        self.calculate_descriptors()


if __name__ == "__main__":
    
    ligand_excel_file = 'ligands_test.xlsx'
    auxiliary_ligands = []
    substrate = ['CC#[N:1]']

    geom = 'SP'
    central_atom = '[Ir+3]'
    names_of_xyz_key = 'Cas'


    # MACE input 
    mace_input = {'bidentate_ligands': ligand_excel_file, 
                  'auxiliary_ligands': auxiliary_ligands, 
                  'names_of_xyz': names_of_xyz_key, 
                  'central_atom': central_atom, 
                  'geom': geom, 
                  'substrate': substrate}

    # ChemXpaX input
    
    current_directory = os.path.join(os.getcwd(), "chemspax")  # use path.join()
    # os.chdir(current_directory)
    path_to_substituents = os.path.join(current_directory, "substituents_xyz") 
    path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")
    substituent_df = pd.read_excel('walphos_chemspax.xlsx').dropna()
    substituent_list = np.array(substituent_df[['R1', 'R2', 'R3', 'R4']])
    print(substituent_list)
    names = substituent_df['Name']

    skeleton_list = list(substituent_df['skeleton'])
    path_to_hand_drawn_skeletons = os.path.join(current_directory, "skeletons")
    path_to_output = os.path.join(current_directory, "complexes")

    chemspax_input = {'skeleton_list' : skeleton_list, 
                    'substituent_list' : substituent_list, 
                    'path_to_database' : path_to_database, 
                    'path_to_substituents' : path_to_substituents, 
                    'path_to_additional_skeletons' : path_to_hand_drawn_skeletons}


    method = 'gfn2'
    charge_of_complex = 0
    multiplicity = 2
    solvent = 'ch2cl2'
    
    crest_input = {'method': method, 
                   'charge_complex': charge_of_complex,
                   'multiplicity' : multiplicity,
                   'solvent' : solvent,
                   'conformer_search' : 'off'}

    # print(skeleton_list)
    workflow = Workflow(chemspax_input = chemspax_input, path_to_workflow = os.getcwd() + '/Workflow')
    
    workflow.run_chemspax(names=names ,skeleton_list_=skeleton_list)
    
