# OBeLiX
![logo](../images/logo.png)
Open Bidentate Ligand eXplorer (OBeLiX) is a workflow for automated and reproducible TM-based catalyst structure generation and descriptor calculation. It Uses
our in-house developed packages Mace and ChemSpaX for
bias-free catalyst structure generation. CREST is used for conformer search and Morfeus/cclib
are used for descriptor extraction. The workflow is currently usable for monodentate and bidentate ligands.

## Installation

The workflow can be installed from the environment.yml file with conda:

```bash
conda env create -f environment.yml
```

The pip installable version will be shortly available.

## Folder structure

The output folder structure is a crucial first step in using the workflow. You can either create your own workflow structure and then run the workflow with your own paths. Or it is recommended that the below commands are used:

```python
workflow = Workflow(path_to_workflow = "your/path/to/workflow_results")
workflow.prepare_folder_structure()
```
The workflow contains 4 modules: MACE, ChemSpaX, xTB/CREST, descriptor calculation.
This command will thus create a folder structure like this at the path you specified:
```bash
Workflow/MACE
Workflow/ChemSpaX
Workflow/CREST
Workflow/Descriptors
```

## Workflow Structure
### Using all modules at once
You can run the whole workflow at once, like this:

```python
workflow = Workflow(mace_input = mace_input, chemspax_input = chemspax_input, crest_input = crest_input, path_to_workflow = "your/path/to/workflow_results")
workflow.run_workflow()
```

The input to the workflow class contains the input to MACE, input to ChemSpax, input to xTB/CREST and the descriptor calculator. The inputs should be given as below:

```python
from run_workflow import *

ligand_excel_file = os.path.join(os.getcwd(), 'test_mace.xlsx')
auxiliary_ligands = ['CC#[N:1]', 'CC#[N:1]']
# substrate = ['CC#[N:1]'] ## for octahedral complexes one can use alternative ligand names.

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
skeleton_list = substituent_df['Functionalization']
# print(skeleton_list)
path_to_hand_drawn_skeletons = os.path.join(current_directory, "skeletons")
path_to_output = os.path.join(current_directory, "complexes")

chemspax_input = {'substituent_list' : substituent_list, 
                'path_to_database' : path_to_database, 
                'path_to_substituents' : path_to_substituents}

workflow = Workflow(mace_input = mace_input, chemspax_input=chemspax_input, path_to_workflow = os.getcwd() + '/wf_test5', geom='BD')
workflow.prepare_folder_structure()
workflow.run_mace()
workflow.run_chemspax(names=names, functionalization_list=skeleton_list)

```

### Using the modules separately


## Citation
WIP