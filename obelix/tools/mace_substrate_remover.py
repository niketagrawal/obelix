from obelix.run_workflow import *
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from morfeus.io import write_xyz


def remove_substrate(substrate_SMILES, ligand_keys,data_path):
    """ This function removes the substrate from MACE generated xyz files.
    args:
    substrate_SMILES: SMILES of substrate in list format
    ligand_keys: name of ligands in list format, can be multiple ligands
    data_path: path to xyz files, can contain multiple xyz files
    metal_center: metal center of complexes in string format, used for convenient naming of new xyz file
    return: xyz files without substrate, denoted wuth the suffix "no_substrate", saved in a new directory called substrate_removed"""

    # Create directory for xyz files without substrate
    if not os.path.exists(os.path.join(data_path, 'substrate_removed')):
        os.makedirs(os.path.join(data_path, 'substrate_removed'))
    
    # Convert substrate SMILES to molecular formula and extract atoms and counts
    substrate_mol = Chem.MolFromSmiles(substrate_SMILES[0])
    substrate_mol = CalcMolFormula(substrate_mol)
    print(f'Molecular formula of substrate: {substrate_mol}')
    substrate_df = pd.DataFrame(columns=['Atom', 'Count'])
    for i in range(len(substrate_mol)):
        if substrate_mol[i].isupper():
            atom = substrate_mol[i]
            count = 1  # Standard count if no digit is found after the atom

            # Check if the next character is a digit, if so, that will be the count of that specific atom
            if i+1 < len(substrate_mol) and substrate_mol[i+1].isdigit():
                count = int(substrate_mol[i+1])
            if i+2 < len(substrate_mol) and substrate_mol[i+1].isdigit() and substrate_mol[i+2].isdigit():
                count = int(substrate_mol[i+1] + substrate_mol[i+2])
                
            # Add new atom and count to substrate dataframe
            substrate_df = pd.concat([substrate_df, pd.DataFrame({'Atom': [atom], 'Count': [count]})], ignore_index=True)
        # If atom consits of multiple characters (e.g Cl), fix atom name and count
        if substrate_mol[i].islower():
            substrate_df['Atom'].iloc[-1] = substrate_mol[i-1]+substrate_mol[i]
            if i+1 < len(substrate_mol) and substrate_mol[i+1].isdigit():
                substrate_df['Count'].iloc[-1] = int(substrate_mol[i+1])
            if i+2 < len(substrate_mol) and substrate_mol[i+1].isdigit() and substrate_mol[i+2].isdigit():
                substrate_df['Count'].iloc[-1] = int(substrate_mol[i+1] + substrate_mol[i+2])
    print(f'Decomposed substrate: \n{substrate_df}')
    # Find number of hydrogen and non-hydrogen atoms in substrate
    H_count = substrate_df['Count'][substrate_df['Atom'] == 'H'].sum()
    non_H_count = sum(substrate_df['Count'])-H_count
    # Import xyz files, remove substrate and write new xyz files
    for ligand in ligand_keys:
        xyz_file = glob.glob(os.path.join(glob.escape(data_path),'*'+ligand+'*.xyz'))
        for f in range(len(xyz_file)):
            x,y = read_xyz(xyz_file[f])

            # Create dataframe from xyz file
            df = pd.DataFrame({'Atom':x, 'x-coordinates':y.transpose()[0], 'y-coordinates':y.transpose()[1], 'z-coordinates':y.transpose()[2]})
            # Find the index where hydrogen atoms connected to carbon start
            for i, row in df.iterrows():
                if row['Atom'] != 'H':
                    hydrogen_index = i+1 
            # Remove substrate H atoms from xyz
            df = pd.concat([df[:hydrogen_index],df[hydrogen_index+H_count:]], ignore_index=True)
            # Remove substrate non-H atoms from xyz
            df = df[non_H_count:]
            df.reset_index()
            # Write new xyz file
            xyz_filename = os.path.basename(xyz_file[f]).replace('.xyz','_no_substrate.xyz')
            write_xyz(os.path.join(data_path, 'substrate_removed', xyz_filename), df['Atom'], df[['x-coordinates', 'y-coordinates', 'z-coordinates']].values.tolist())
    return print(f"Substrate removed for {len(ligand_keys)} ligands")

if __name__ == '__main__':
    os.chdir(os.path.normpath(os.getcwd() + os.sep + os.pardir))
    data_path = os.path.join(os.getcwd(),'acetonitrile', '[Ru+2]', 'PP', 'MACE','H-H_axial')

    smiles_path = os.path.join(os.getcwd(),'ligand_smiles','ligand_smiles_non_ferrocenes.xlsx')
    smiles_df = pd.read_excel(smiles_path).dropna()
    ligand_number = smiles_df['Ligand#']

    remove_substrate(['CC#[N:1]'], ligand_number, data_path)
