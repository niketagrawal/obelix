from chemspax.main import main
from chemspax.utilities import *
import os
import glob
import pandas as pd
from tools.utilities import *s
from tools.chemspax_util import *
import shutil
from tqdm import tqdm
import mace
import subprocess
import morfeus as mf
from morfeus.conformer import ConformerEnsemble
from morfeus import BiteAngle, ConeAngle, BuriedVolume, Dispersion, SASA, read_xyz
from xtb.utils import get_method
import pandas as pd
from tools.chemspax_util import prepare_data
from openbabel import openbabel


class MACE:   
    def __init__(self, bidentate, CA, name_of_xyz):
        # Declare instances of the class (a.k.a self variable)
        self.CA = CA
        self.name_of_xyz = name_of_xyz
        self.bidentate = bidentate
    
    def generate_complex_SP_xyz(self):
        geom = 'SP'

        core = mace.ComplexFromLigands([self.bidentate], self.CA, geom)
        Xs = core.GetStereomers(regime='all', dropEnantiomers=True)
        for i, X in enumerate(Xs):
            X.AddConformers(numConfs=10)   
            X.ToXYZ(self.CA + '_' + '{}{}.xyz'.format(self.name_of_xyz, i), confId='min')
            
        
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
        for i, X in enumerate(Xs):
            X.AddConformers(numConfs=10)   
            X.ToXYZ(self.CA + '_' + '{}{}.xyz'.format(self.name_of_xyz, i), confId='min')
            


class Workflow:
    def __init__(self, mace_input = [], chemspax_input = [], crest_input = [], path_to_workflow = []):
        self.mace_input = mace_input
        self.chemspax_input = chemspax_input
        self.crest_input = crest_input
        self.path_to_workflow = path_to_workflow
        
        
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
        names_of_xyz_key = list(pd.read_excel(self.mace_input['bidentate_ligands'])['Cas'])
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
                  complex.generate_complex_OH_xyz(auxiliary_ligands, substrate=self.substrate) 
              
              if self.geom == 'SP':
                  identifier_SP = self.names_of_xyz[i] + '_SP'
                  complex = MACE(self.mace_ligands[i], central_atom,  identifier_SP)
                  complex.generate_complex_SP_xyz()
            except Exception:
              print('Wrong SMILES:', i)               
    
    def run_chemspax(self):  
        
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
            main(skeleton_list, sub_list, self.path_to_database, path_to_substituents, path_to_skeletons, chemspax_working_directory, path_to_output)
            os.rename(os.path.join(path_to_output, os.path.basename(os.path.normpath(skeleton_list[0]))[:-4] + '_func_' + str(len(sub_list)) + '.mol'), \
                os.path.join(path_to_output, os.path.basename(os.path.normpath(skeleton_list[0]))[:-4] + '_' + '_'.join(sub_list) + '.mol'))
            
            ## Convert mol to xyz to keep the bonding information
            obconversion = openbabel.OBConversion()
            obconversion.SetInFormat('mol')
            obconversion.SetOutFormat('xyz')
            mol = openbabel.OBMol()
            obconversion.ReadFile(mol, \
                os.path.join(path_to_output, os.path.basename(os.path.normpath(skeleton_list[0]))[:-4] + '_' + '_'.join(sub_list) + '.mol'))
            
            obconversion.WriteFile(mol, \
                os.path.join(path_to_output, os.path.basename(os.path.normpath(skeleton_list[0]))[:-4] + '_' + '_'.join(sub_list) + '.xyz'))
            
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
            
            try:
                for complex in complexes_to_calc_descriptors:
                    properties = {}
                    elements, coordinates = read_xyz(os.path.join(complex, 'xtbopt.xyz'))
                    
                    
                    bidentate = find_bidentate(os.path.join(complex, 'xtbopt.xyz'))
                    properties["bite_angle"] = BiteAngle(coordinates, bidentate[1], bidentate[0], bidentate[2]).angle
                    
                    # elements_no_H = []
                    # H_indices = []
                    # for elem_index, elem in enumerate(elements):
                    #     if elem != 'H':
                    #         elements_no_H.append(elem)     
                    #     else:
                    #         H_indices.append(elem_index)
                    # H_indices.sort(reverse=True)
                    # coordinates_no_H = np.delete(np.array(coordinates), H_indices, axis=0)
                    
                    # nr_of_atoms_less_than_metal = 0

                    # for index in H_indices:
                    #     if index > bidentate[1]:
                    #         pass
                    #     else:
                    #         nr_of_atoms_less_than_metal += 1
                    
                    
                    properties["cone_angle"] = ConeAngle(elements, coordinates, bidentate[1]).cone_angle
                    
                    # print(properties["cone_angle"])
                    bv1 = BuriedVolume(elements, coordinates, bidentate[1], radius=3.5)
                    bv2 = BuriedVolume(elements, coordinates, bidentate[0], radius=3.5)
                    bv3 = BuriedVolume(elements, coordinates, bidentate[2], radius=3.5)
                    
                    properties["buried_volume_metal_center"] = bv1.fraction_buried_volume
                    properties["buried_volume_P1"] = bv2.fraction_buried_volume
                    properties["buried_volume_P2"] = bv3.fraction_buried_volume

                    # print(BuriedVolume(ce.elements, conformer.coordinates, bidentate[1]).print_report())
                    # BuriedVolume(elements, coordinates, bidentate[1], radius=4).print_report()

                    properties["dispersion"] = Dispersion(elements, coordinates).atom_p_int[bidentate[1]]
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
                    properties["HOMO_LUMO_gap"] = homo - lumo
                    
                    for property in properties.keys():
                            dictionary_for_properties[os.path.basename(os.path.normpath(complex))] = properties
            except Exception:
                print("Descriptor calculation failed for this complex:", os.path.basename(os.path.normpath(complex)))   
                         
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
                    
                    bidentate = find_bidentate(conformer)
                    conformer.properties["bite_angle"] = BiteAngle(conformer.coordinates, bidentate[1], bidentate[0], bidentate[2]).angle
                    H_indices = []
                    
                    for elem_index, elem in enumerate(elements):
                        if elem == 'H':
                            H_indices.append(elem_index)                    

                    coordinates_no_H = np.delete(np.array(coordinates), H_indices, axis=0)
                    
                    elements_no_H = elements[:]
                    
                    for elem_index in H_indices:
                        H_indices.pop(elem_index)
                    
                    conformer.properties["cone_angle"] = ConeAngle(elements_no_H, coordinates_no_H, bidentate[1]).cone_angle
                    conformer.properties["buried_volume"] = BuriedVolume(ce.elements, conformer.coordinates, bidentate[1]).fraction_buried_volume
                    conformer.properties["dispersion"] = Dispersion(ce.elements, conformer.coordinates).atom_p_int[1]
                    conformer.properties["sasa"] = SASA(ce.elements, conformer.coordinates).area
                    
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
        dataframe.to_excel('descriptors_with_specific_buried_vol.xlsx')     
        
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
    substituent_list = [["C6H12", "C6H12", "3_5_trifluoro_methyl_phenyl", "3_5_trifluoro_methyl_phenyl"]]
   
    # substituent_list = [["CCH3CH3CH3", "CCH3CH3CH3", "para_trifluoromethyl_phenyl", "para_trifluoromethyl_phenyl"], 
    #                     ["C6H12", "C6H12", "para_trifluoromethyl_phenyl", "para_trifluoromethyl_phenyl"], 
    #                     ["C6H12", "C6H12", "3_5_dimethyl_4_metoxy_phenyl", "3_5_dimethyl_4_metoxy_phenyl"],
    #                     ["CCH3CH3CH3", "CCH3CH3CH3", "furanyl", "furanyl"],
    #                     ["3_5_dimethyl_phenyl", "3_5_dimethyl_phenyl", "3_5_dimethyl_4_metoxy_phenyl", "3_5_dimethyl_4_metoxy_phenyl"],
    #                     ["C6H6-CH3_ortho_1", "C6H6-CH3_ortho_1", "furanyl", "furanyl"],
    #                     ["C6H6-CH3_ortho_1", "C6H6-CH3_ortho_1", "CCH3CH3CH3", "CCH3CH3CH3"], 
    #                     ["3_5_dimethyl_phenyl", "3_5_dimethyl_phenyl", "3_5_trifluoro_methyl_phenyl", "3_5_trifluoro_methyl_phenyl"],
    #                     ["C6H6", "C6H6", "C6H12", "C6H12"],
    #                     ["3_5_dimethyl_phenyl", "3_5_dimethyl_phenyl", "C6H6", "C6H6"], 
    #                     ["CCH3CH3CH3", "CCH3CH3CH3", "3_5_dimethyl_4_metoxy_phenyl", "3_5_dimethyl_4_metoxy_phenyl"],
    #                     ["3_5_dimethyl_phenyl", "3_5_dimethyl_phenyl", "naphtalenyl", "naphtalenyl"]]  # will find substituent.xyz
    

    # substituent_list = [["3_5_dimethyl_phenyl", "3_5_dimethyl_phenyl", "3_5_trifluoro_methyl_phenyl", "3_5_trifluoro_methyl_phenyl"]]
    skeleton_list = glob.glob(os.path.join("skeletons", "*.xyz"))
    path_to_hand_drawn_skeletons = os.path.join(current_directory, "skeletons")
    path_to_output = os.path.join(current_directory, "complexes")


    chemspax_input = {'skeleton_list' : skeleton_list, 
                    'substituent_list' : substituent_list, 
                    'path_to_database' : path_to_database, 
                    'path_to_substituents' : path_to_substituents, 
                    'path_to_additional_skeletons' : path_to_hand_drawn_skeletons}


    method = 'gfn2'
    charge_of_complex = 0
    multiplicity = 1
    solvent = 'ch2cl2'
    
    crest_input = {'method': method, 
                   'charge_complex': charge_of_complex,
                   'multiplicity' : multiplicity,
                   'solvent' : solvent,
                   'conformer_search' : 'off'}

    workflow = Workflow(mace_input = mace_input, chemspax_input = chemspax_input, path_to_workflow = os.getcwd() + '/Workflow')
    workflow.calculate_descriptors(path_to_outputs = os.path.join(os.getcwd(), "Workflow", "CREST"), output = 'xtb')
    
    
