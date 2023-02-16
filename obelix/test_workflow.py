from run_workflow import *

ligand_excel_file = os.path.join(os.getcwd(), 'test_mace.xlsx')
auxiliary_ligands = ['CC#[N:1]', 'CC#[N:1]']
substrate = ['CC#[N:1]']

geom = 'SP'
central_atom = '[Rh+]'
names_of_xyz_key = 'Cas'


# MACE input 
mace_input = {'bidentate_ligands': ligand_excel_file, 
                'auxiliary_ligands': auxiliary_ligands, 
                'names_of_xyz': names_of_xyz_key, 
                'central_atom': central_atom, 
                'geom': geom, 
                'substrate': substrate}

current_directory = os.path.join(os.getcwd(), "chemspax")  # use path.join()
# os.chdir(current_directory)
path_to_substituents = os.path.join(current_directory, "substituents_xyz") 
path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")
substituent_df = pd.read_excel('test_chemspax.xlsx').dropna()
substituent_list = np.array(substituent_df[['R1', 'R2', 'R3', 'R4']])
# print(substituent_list)
names = substituent_df['Name']
func_list = substituent_df['Functionalization']
# print(skeleton_list)
path_to_hand_drawn_skeletons = os.path.join(current_directory, "skeletons")
path_to_output = os.path.join(current_directory, "complexes")

chemspax_input = {'substituent_list' : substituent_list, 
                'path_to_database' : path_to_database, 
                'path_to_substituents' : path_to_substituents}

workflow = Workflow(mace_input = mace_input, chemspax_input=chemspax_input, path_to_workflow = os.getcwd() + '/wf_test5', geom='BD')
workflow.prepare_folder_structure()
workflow.run_mace()
workflow.run_chemspax(names=names, functionalization_list=func_list)
