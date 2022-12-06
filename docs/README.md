# open-bidentate-explorer

A workflow for automated catalyst structure analysis. Uses Mace and ChemSpaX for
bias-free catalyst structure generation. CREST is used for conformer search and Mofeus/Mordred
are used for descriptor calculation. Afterwards, the descriptors are boltzmann averaged and compared for
predictive purposes.


# Folder structure

The folder structure is a crucial first step in using the workflow. You can either create your own workflow structure and then run the workflow with your own paths. Or alternatively you could use:

```
workflow = Workflow(path_to_workflow = 'your/path/to/workflow")
workflow.prepare_folder_structure()
```
