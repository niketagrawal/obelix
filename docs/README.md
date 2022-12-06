# open-bidentate-explorer

A workflow for automated catalyst structure analysis. Uses Mace and ChemSpaX for
bias-free catalyst structure generation. CREST is used for conformer search and Mofeus/Mordred
are used for descriptor calculation. Afterwards, the descriptors are boltzmann averaged and compared for
predictive purposes.

## Installation

The workflow can be installed from the environment.yml file with conda:

```bash
conda env create -f environment.yml
```

The pip installable version will be shortly available.

## Folder structure

The folder structure is a crucial first step in using the workflow. You can either create your own workflow structure and then run the workflow with your own paths. Or alternatively you could use:

```python
workflow = Workflow(path_to_workflow = "your/path/to/workflow")
workflow.prepare_folder_structure()
```

## Workflow Structure

The workflow contains 4 modules: MACE, ChemSpaX, xTB/CREST, Morfeus.

You can run the whole workflow at once, like this: 

```python
workflow = Workflow(mace_input = mace_input, chemspax_input = chemspax_input, crest_input = crest_input, path_to_workflow = "your/path/to/workflow")
workflow.prepare_folder_structure()
```


The input to the workflow class contains the input to MACE, input to ChemSpax, input to xTB/CREST. The inputs should be given as below:

```python
ligand_excel_file = 'ligands_test.xlsx'
auxiliary_ligands = []
substrate = ['CC#[N:1]']

geom = 'SP'
central_atom = '[Ir+3]'
names_of_xyz_key = 'Cas' # Names attributed to the MACE generated files


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
workflow.run_workflow()

```
