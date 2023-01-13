import random, string
from rdkit import Chem
import numpy as np
import pandas as pd
import glob, os


def add_code_to_structure():
    N = 10
    code_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))
    return code_name
    

def change_second_line_xyz(xyz, new_content=''):
    with open('{}'.format(xyz), 'r', newline='\n') as file:
        # read a list of lines into data
        data = file.readlines()
    # now change the 2nd line, note that you have to add a newline
    
    if new_content == '':
        data[1] = '\n'
    else:
        data[1] = '{}\n'.format(new_content)
        
    # and write everything back
    with open('{}'.format(xyz), 'w') as file:
        file.writelines(data)


def get_ligands_from_smiles(ligands_SMILES):
    ligand_mol = ['']*len(ligands_SMILES)
    for index, ligand in enumerate(ligands_SMILES):
        ligands_SMILES[index] = Chem.MolFromSmiles(ligand)
    
    return ligand_mol


def dataframe_from_dictionary(dictionary):
    Dataframe = pd.DataFrame.from_dict({i: dictionary[i]
                           for i in dictionary.keys()}, orient='index')
    return Dataframe


def get_cluster_centroid_coord(n_clusters, unique_labels, dataframe):

    """
    This function identifies the centroids of the clusters generated by the UMAP 
    algorithm. 
    
    Args:
        n_clusters    : integer
        unique_labels : array of integers 
        dataframe     : pandas dataframe containing UMAP1, UMAP2
    Return:
        centroids     : array of float sized (2, n_clusters), 
                          where 2 is (x, y) coordinates
    """

    # initialize cluster centroids array
    centroids = np.zeros((2, n_clusters))
    
    # fill the array with the coords
    for j, u_label in enumerate(unique_labels):
        UMAP1 = dataframe[(dataframe['label'] == u_label)].mean()['UMAP-1']
        UMAP2 =  dataframe[(dataframe['label'] == u_label)].mean()['UMAP-2']    
        temp_coord_array = np.array([UMAP1, UMAP2])
        centroids[1][j] = temp_coord_array[0]
        centroids[0][j] = temp_coord_array[1]
        
    return centroids


def find_bonds_with_neighbours(mol, center_atom_atomic_number):
    center_atom = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == center_atom_atomic_number][0]
    bonds = [bond for bond in center_atom.GetBonds()]
    return bonds


def check_if_at_least_two_mapped_atoms_in_ring(list_of_mapped_idxs, list_of_ring_idxs):
    return len(set(list_of_mapped_idxs) & set(list_of_ring_idxs)) >= 2


def find_bidentate(xyz):

    with open(xyz) as file:
        listfile = []

        for line in file:
            listfile.append(line.strip())
    indices = [0]*3
    atoms = []
    for index, atom in enumerate(listfile):
        if atom[0] == 'P':
            if not indices[0]:
                indices[0] = index - 1
            else:
                indices[2] = index - 1
            atoms.append(atom[0])
        else:
            if atom[:2] == 'Ir':
                indices[1] = index - 1
                atoms.append(atom[:2])
                
    return indices


def xyz_to_gjf(header, io_path):
    
    if not os.path.exists(io_path):
        os.mkdir(io_path)

    os.chdir(io_path)

    xyzs = glob.glob('*.xyz')
    
    for xyz in xyzs:
        with open(xyz) as file:
            f = file.readlines()
            print(len(f))
            for new_line in reversed(header):
                f.insert(0, new_line)        
            
            f.pop(len(header))
            f.pop(len(header))
            
            if f[-1] != '\n':
                f.append('\n')
            file.close()
  
            with open(f'{xyz[:-4]}.gjf', 'w+') as file:
            
                file.writelines(f)
                file.close()


def gjf_to_xyz(path, header):
    for gjf in glob.glob(os.path.join(path, "*.gjf")):
        with open(gjf) as file:
            f = file.readlines()
            
            # count = len(f) - len(header)
            # print(count, gjf)
            
            f = f[len(header):]
            # print(f)
            count = len(f) - 2
            f.insert(0, '\n')
            f.insert(0, str(count))
            file.close()
        with open(f'{gjf[:-4]}.xyz', 'w+') as file:
            file.writelines(f)
            file.close()  


if __name__ == "__main__":
    header = ['%nprocshared=32 \n',  
            '%mem=32GB \n', 
            '#p opt freq pop=nbo def2svpp empiricaldispersion=gd3bj integral=grid=ultrafinegrid pbe1pbe\n', 
            'scf=(xqc,maxconventionalcycles=90) nosymm\n',      
            '\n', 
            'Title Card Required\n'
            '\n', 
            '1 1 \n']


    # xyz_to_gjf(header=header, 
    #            io_path='data/Rh_dataset/ferrocenes/other')


    # gjf_to_xyz(path='data/Rh_dataset/ferrocenes/gaussian_input_files/gjf', header=header)
    xyz_to_gjf(header=header, io_path='../')