import numpy as np 
from numpy.linalg import norm as cartesian_distance
from morfeus import read_xyz
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Needed data

atom_covalent_bond_nr = {'N'  : 2, 
                         'P'  : 2, 
                         'S'  : 2.2, 
                         'C'  : 2.1, 
                         'O'  : 1.8,
                         'H'  : 1.15,
                         'Fe' : 2.1}

metal_centers = {'Rh' : {'OH' : 6, 'BD' : 2}, 
                 'Ir' : {'OH' : 6, 'BD' : 2}, 
                 'Mn' : {'OH' : 6, 'BD' : 2}, 
                 'Rh' : {'OH' : 6, 'BD' : 2},
                 'Ru' : {'OH' : 5, 'BD' : 2}}


def bfs(visited, graph, node):
    queue = []
    
    visited.append(node)
    queue.append(node)

    while queue:
        s = queue.pop(0)
        for neighbour in graph[s]:
            if neighbour not in visited:
                visited.append(neighbour)
                queue.append(neighbour)
    return visited


def molecular_graph(xyz, geom = 'OH'):
    
    # Read xyz - coord, type with morfeus
    elements, coords = read_xyz(xyz)
    
    nr_of_atoms = len(elements)
    
    # init dictionary 
    elem_dict = {}
    
    # Populate dictionary -> index of atom  : {atom type : atom coords}
    
    for index, elem in enumerate(elements):
        elem_dict[index] = {elem: coords[index]}

    # Find atom type like this: atom_type = elem_dict[atom_index].keys() 
    
    # atom_type = elem_dict[16].keys()
    # test_bond = atom_covalent_bond_nr[list(atom_type)[0]]

    interatomic_distances = {}
    
    for atom_key1 in range(0, nr_of_atoms):
        interatomic_distances[atom_key1] = []
        for atom_key2 in range(0, nr_of_atoms):
            atom_type1 = list(elem_dict[atom_key1].keys())[0]
            atom_type2 = list(elem_dict[atom_key2].keys())[0]
            interatomic_distances[atom_key1].append(cartesian_distance(elem_dict[atom_key1][atom_type1] - elem_dict[atom_key2][atom_type2]))

    for atom_key in range(0, nr_of_atoms + 1):                

        atom_type = list(elem_dict[atom_key].keys())[0]

        if atom_type in list(metal_centers.keys()):
            nr_of_bonds = metal_centers[atom_type][geom]

            # Arange shortest bonds
            idx = np.argpartition(interatomic_distances[atom_key], nr_of_bonds + 1)
            idx = idx[:nr_of_bonds + 1]
            metal_key = atom_key
            ligand_start_idx = np.setdiff1d(idx, np.array([atom_key, metal_key]))
            
            # ligand_connected_atoms = zip(ligand_start_idx, elements[ligand_start_idx])
            # print(list(ligand_connected_atoms))
            # print(ligand_start_idx)
            break
    
    store_atoms = []
    store_atoms.extend(ligand_start_idx)
    mol_graph = {}

    for (atom_ligand, atom_type) in zip(ligand_start_idx, elements[ligand_start_idx]):      
        
        if atom_type in list(atom_covalent_bond_nr.keys()):
            nr_of_bonds = atom_covalent_bond_nr[atom_type] 
            # print(atom_type)
            idx = np.where(np.array(interatomic_distances[atom_ligand]) <= atom_covalent_bond_nr[atom_type])
            # print('idx', idx)
            idx = idx[:len(idx) + 1]

            # bond_distances = np.array(interatomic_distances[atom_ligand])[idx]                           
            bonded_atoms = list(np.setdiff1d(idx, np.array([atom_ligand, metal_key])))
            # bonded_atoms.append(metal_key)
            store_atoms.extend(bonded_atoms)
            mol_graph[atom_ligand] = bonded_atoms
    
    ligandMetal = list(ligand_start_idx)
    ligandMetal.append(metal_key)
    # print('Ligand_metal', ligandMetal)
    
    for atom_key in range(0, nr_of_atoms):
        if atom_key not in ligandMetal:
            atom_type = elements[atom_key]
            nr_of_bonds = atom_covalent_bond_nr[atom_type] 
             
            idx = np.where(np.array(interatomic_distances[atom_key]) <= atom_covalent_bond_nr[atom_type])
            idx = idx[:len(idx) + 1]

            # bond_distances = np.array(interatomic_distances[atom_ligand])[idx]                           
            bonded_atoms = list(np.setdiff1d(idx, np.array([atom_key])))
            
            # for id_atom, atom in enumerate(bonded_atoms):
            #     if atom in store_atoms:
            #         del bonded_atoms[id_atom]
                    
            store_atoms.extend(bonded_atoms)
            mol_graph[atom_key] = bonded_atoms
            
    mol_graph[metal_key] = ligand_start_idx

    ### Adjacency Matrix (fingerprint)
    # keys=sorted(mol_graph.keys())
    # size=len(keys)

    # adj_matrix = [ [0]*size for i in range(size) ]

    # for a,b in [(keys.index(a), keys.index(b)) for a, row in mol_graph.items() for b in row]:
    #     adj_matrix[a][b] = 1
        
    ### Breadth first search algorithm (go through all neighbours and store)
    ligands_atoms_idx = {}
    checklist = []
    for atom_ in ligand_start_idx:
        visited = []
        ligands_atoms_idx[atom_] = bfs(visited, mol_graph, atom_)    
        checklist.extend(ligands_atoms_idx[atom_])

    # Plus one because the metal is excluded to find the subgraphs
    print(f'Length of checklist: ', {len(np.unique(np.array(checklist))) + 1}, 'Nr. of atoms: ', {len(elements)})
    # print(np.array(adj_matrix).shape)
    # print(np.unique(np.array(store_atoms).sort))
    # To find bidentate see if any of the ligands map into each other
    
    bidentate_indices = [metal_key]
    
    for ligand_index1 in ligand_start_idx:
        for ligand_index2 in ligand_start_idx:
            if ligand_index1 != ligand_index2 and ligand_index1 in ligands_atoms_idx[ligand_index2]:
                bidentate_indices.append(ligand_index1)
                bidentate_indices.append(ligand_index2)
            break
    
    for elem in ligands_atoms_idx:
        print(elements[elem]) 
        
    return ligands_atoms_idx,bidentate_indices

                    
# print(molecular_graph('xtbopt.xyz'))

ligands_indices, bidentate = molecular_graph('xtbopt.xyz')
print(bidentate, ligands_indices)
