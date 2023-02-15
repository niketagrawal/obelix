from run_workflow import *

ligand_excel_file = os.path.join(os.getcwd(), 'ligand_test.xlsx')
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

current_directory = os.path.join(os.getcwd(), "chemspax")  # use path.join()
# os.chdir(current_directory)
path_to_substituents = os.path.join(current_directory, "substituents_xyz") 
path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")
substituent_df = pd.read_excel('demo_chemspax_cv.xlsx').dropna()
substituent_list = np.array(substituent_df[['R1', 'R2', 'R3', 'R4']])
# print(substituent_list)
names = substituent_df['Name']
skeleton_list = substituent_df['Skeleton']
print(skeleton_list)
# print(skeleton_list)
path_to_hand_drawn_skeletons = os.path.join(current_directory, "skeletons")
path_to_output = os.path.join(current_directory, "complexes")

chemspax_input = {'skeleton_list' : skeleton_list, 
                'substituent_list' : substituent_list, 
                'path_to_database' : path_to_database, 
                'path_to_substituents' : path_to_substituents, 
                'path_to_additional_skeletons' : path_to_hand_drawn_skeletons}

workflow = Workflow(mace_input = mace_input, chemspax_input=chemspax_input, path_to_workflow = os.getcwd() + '/wf_test2', geom='BD')
workflow.prepare_folder_structure()
workflow.run_mace()
workflow.run_chemspax(names=names, skeleton_list_=skeleton_list)
