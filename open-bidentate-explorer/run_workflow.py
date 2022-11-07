from chemspax.main import main
from chemspax.utilities import *
import os
import glob
import pandas as pd
from tools.utilities import *
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
        print(auxiliary)
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
        print('Workflow is initializing. Converting your graph input to variables.')
        print('')
        
        if mace_input != []:  
          self.mace_ligands, self.auxiliary_ligands, self.geom, self.central_atom, self.names_of_xyz, self.substrate = self.initialize_mace()
        if chemspax_input != []:        
          self.substituent_list, self.path_to_database, self.path_to_substituents = self.initialize_chemspax()
        if crest_input != []:
          self.method, self.charge_of_complex, self.multiplicity, self.solvent = self.initialize_crest()
                  
          

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
        
        return method, charge_of_complex, multiplicity, solvent
        
        
    def prepare_folder_structure(self):
        
        print('Preparing the folder structure')
    
        os.mkdir(self.path_to_workflow)
        os.chdir(self.path_to_workflow)
        os.mkdir('MACE')
        os.mkdir('ChemSpaX')
        os.mkdir('CREST')
        os.mkdir('Descriptors')      
        os.chdir('..')


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
        os.chdir("..")
        
        
    def run_chemspax(self, functionalization_list = ''):  
        
        """
        Run chemspax from the pip installable package chemspax.
        Use prepare_data(...) to prepare the substituents.
        Then run main(...) to functionalize the skeletons prepared by MACE.
        
        :param functionalization_list -> this will go to the second line of the copied MACE skeletons into the ChemSpax/skeletons dir.
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
          skeleton_list = glob.glob(self.path_to_workflow + '/ChemSpaX/skeletons/*.xyz')
        
          for xyz_file in skeleton_list:
              change_second_line_xyz(xyz_file, new_content=functionalization_list)
  
          # set skeleton path to ../ChemSpaX/skeletons        
          path_to_skeletons = os.path.join(os.getcwd(), dest_dir_skeletons)
          
          # copy substituents from user source to ChemSpaX folder
          src_dir_subs = self.path_to_substituents
          dest_dir_subs = os.path.join(self.path_to_workflow, 'ChemSpaX', 'substituents_xyz')        
          shutil.copytree(src_dir_subs, dest_dir_subs)
  
          # Prepare data (from data_preparation.py in the chemspax package)
                  
          print('Data preparation has been performed.')
          path_to_substituents = os.path.join(os.getcwd(), dest_dir_subs)
          
          # Set path to output
          path_to_output = os.path.join(self.path_to_workflow, 'ChemSpaX','chemspax_output')
                  
          # Run chemspax from main
        for index, sub_list in enumerate(self.substituent_list):
            main(skeleton_list, sub_list, path_to_database, path_to_substituents, path_to_skeletons, chemspax_working_directory, path_to_output)
            os.rename(os.path.join(path_to_output, os.path.basename(os.path.normpath(skeleton_list[0]))[:-4] + '_func_' + str(len(sub_list)) + '.xyz'), \
                os.path.join(path_to_output, os.path.basename(os.path.normpath(skeleton_list[0]))[:-4] + '_' + '_'.join(sub_list) + '.xyz'))
            for filename in glob.glob(os.path.join(self.path_to_workflow, 'ChemSpaX', 'chemspax_output', '*_func_*')):
                os.remove(filename) 
    
    def run_crest(self, path_to_complexes = []):
        
        """
        The crest call takes 6 inputs:
        
        path_to_complexes=$1
        method=$2
        charge=$3
        geom=$4
        multiplicity=$5
        solvent=$6
        """
    
        print('About to run optimization:')
        
        if path_to_complexes == []:
          path_to_complexes = os.path.join(self.path_to_workflow, 'ChemSpaX', 'chemspax_output')
          os.chdir(self.path_to_workflow)
          
        subprocess.call('bash "../run_crest.sh" "%s" %s %s %s %s %s' % (path_to_complexes, self.method, self.charge_of_complex, self.geom, self.multiplicity, self.solvent), shell=True)
        os.chdir("-")
    
    def calculate_descriptors(self, path_to_outputs = [], output = 'CREST'):
        
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
                elements, coordinates = read_xyz(os.path.join(complex, 'xtbopt.xyz'))
                
                bidentate = find_bidentate(os.path.join(complex, 'xtbopt.xyz'))
                properties["bite_angle"] = BiteAngle(coordinates, bidentate[1], bidentate[0], bidentate[2]).angle
                properties["cone_angle"] = ConeAngle(elements, coordinates, bidentate[1]).cone_angle
                properties["buried_volume"] = BuriedVolume(elements, coordinates, bidentate[1]).fraction_buried_volume
                
                properties["dispersion"] = Dispersion(elements, coordinates).atom_p_int[1]
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
                    conformer.properties["cone_angle"] = ConeAngle(ce.elements, conformer.coordinates, bidentate[1]).cone_angle
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
        dataframe.to_excel('descriptors_NEW.xlsx')     
        
           
    def run_workflow(self):
        
        if self.path_to_workflow != []:
           self.prepare_folder_structure()
        
        self.run_mace()
        self.run_chemspax()
        self.run_crest()
        self.calculate_descriptors(output='xtb')


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
    substituent_list = [["CCH3CH3CH3", "CCH3CH3CH3", "para_trifluoromethyl_phenyl", "para_trifluoromethyl_phenyl"], 
                        ["C6H12", "C6H12", "para_trifluoromethyl_phenyl", "para_trifluoromethyl_phenyl"], 
                        ["C6H12", "C6H12", "3_5_dimethyl_4_metoxy_phenyl", "3_5_dimethyl_4_metoxy_phenyl"],
                        ["CCH3CH3CH3", "CCH3CH3CH3", "furanyl", "furanyl"],
                        ["3_5_dimethyl_phenyl", "3_5_dimethyl_phenyl", "3_5_dimethyl_4_metoxy_phenyl", "3_5_dimethyl_4_metoxy_phenyl"],
                        ["ortho_methyl_phenyl", "ortho_methyl_phenyl", "furanyl", "furanyl"],
                        ["ortho_methyl_phenyl", "ortho_methyl_phenyl", "CCH3CH3CH3", "CCH3CH3CH3"], 
                        ["3_5_dimethyl_phenyl", "3_5_dimethyl_phenyl", "3_5_triflluoromethyl_phenyl", "3_5_triflluoromethyl_phenyl"],
                        ["C6H6", "C6H6", "C6H12", "C6H12"],
                        ["3_5_dimethyl_phenyl", "3_5_dimethyl_phenyl", "C6H6", "C6H6"], 
                        ["CCH3CH3CH3", "CCH3CH3CH3", "3_5_dimethyl_4_metoxy_phenyl", "3_5_dimethyl_4_metoxy_phenyl"],
                        ["3_5_dimethyl_phenyl", "3_5_dimethyl_phenyl", "naphtalenyl", "naphtalenyl"]]  # will find substituent.xyz
    skeleton_list = glob.glob(os.path.join("skeletons", "*.xyz"))
    path_to_hand_drawn_skeletons = os.path.join(current_directory, "skeletons")
    path_to_output = os.path.join(current_directory, "complexes")
    print(path_to_hand_drawn_skeletons)

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
                   'solvent' : solvent}

    workflow = Workflow(mace_input = mace_input, crest_input = crest_input, chemspax_input=chemspax_input, path_to_workflow = os.getcwd() + '/Workflow1')
    workflow.run_chemspax()
    # workflow.calculate_descriptors(path_to_outputs=os.getcwd() + '/Workflow1/CREST' , output='xtb')
    
    # workflow.run_crest(path_to_complexes = os.path.join(os.getcwd(), "Workflow_test2_SP", "MACE"))
    # workflow.calculate_descriptors(path_to_outputs = os.path.join(os.getcwd(), "Workflow_test2_SP", "CREST"))

