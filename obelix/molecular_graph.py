import numpy as np 
from numpy.linalg import norm as cartesian_distance
from morfeus import read_xyz
from morfeus.utils import convert_elements
import pandas as pd
import periodictable

# Needed data
periodic_table = [element.symbol for element in periodictable.elements]

atom_covalent_max_dist = {}
for element in periodictable.elements:
    if element.interatomic_distance is not None:
        # to be sure that the distance is not too small, we use the interatomic distance
        atom_covalent_max_dist[element.symbol] = element.interatomic_distance
    else:
        # if the interatomic distance is not available, we use a default value
        atom_covalent_max_dist[element.symbol] = 2.5

metal_centers = {'Rh' : {'OH' : 6, 'BD' : 2}, 
                 'Ir' : {'OH' : 6, 'BD' : 2}, 
                 'Mn' : {'OH' : 6, 'BD' : 2}, 
                 'Ru' : {'OH' : 6, 'BD' : 2},
                 'Pd' : {'OH' : 6, 'BD' : 2}}


donor_atoms = ['P', 'N']


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


def molecular_graph(elements, coords, geom = 'OH'):
    
    # Read xyz - coord, type with morfeus
    # elements, coords = read_xyz(xyz)
    if str(elements[0]).isdigit():
        # convert numerical values to strings
        elements = convert_elements(elements, output='symbols')
        # new_elements_mapping = []
        # for index, element in enumerate(elements):
        #     new_elements_mapping.append(periodic_table[element - 1])
        # elements = np.array(new_elements_mapping)
    elements = np.array(elements)
    nr_of_atoms = len(elements)
    # print(nr_of_atoms)
    # init dictionary 
    elem_dict = {}
    
    # Populate dictionary -> index of atom  : {atom type : atom coords}

    for index, elem in enumerate(elements):
        elem_dict[index] = {elem: coords[index]}
    # Find atom type like this: atom_type = elem_dict[atom_index].keys() 
    
    # atom_type = elem_dict[16].keys()
    # test_bond = atom_covalent_max_dist[list(atom_type)[0]]

    interatomic_distances = {}
    
    for atom_key1 in range(0, nr_of_atoms):
        interatomic_distances[atom_key1] = []
        for atom_key2 in range(0, nr_of_atoms):
            atom_type1 = list(elem_dict[atom_key1].keys())[0]
            atom_type2 = list(elem_dict[atom_key2].keys())[0]
            interatomic_distances[atom_key1].append(cartesian_distance(elem_dict[atom_key1][atom_type1] - elem_dict[atom_key2][atom_type2]))

    for atom_key in range(0, nr_of_atoms + 1):
        atom_type = list(elem_dict[atom_key].keys())
        atom_type = atom_type[0]

        if atom_type in list(metal_centers.keys()):
            nr_of_bonds = metal_centers[atom_type][geom]

            # find indices of shortest bonds of donor elements
            metal_key = int(atom_key)
            bidentate_indices = [metal_key]
            # print('geeeeeom', geom)
            if geom == 'OH':
                idx = np.argpartition(interatomic_distances[atom_key], nr_of_bonds + 1)            
                ligand_start_idx = np.setdiff1d(idx, np.array([atom_key, metal_key]))
            else:
                store_donor_atoms = []
                for ind, el in enumerate(elements):
                    if el in donor_atoms:
                        store_donor_atoms.append(ind)
                if len(store_donor_atoms) == 2: 
                    ligand_start_idx = np.array(store_donor_atoms)
                    bidentate_indices.extend(store_donor_atoms)
                    # print(bidentate_indices)
                    break
                else:
                    # goal is to get the 2 shortest donor atoms.
                    # print(np.array(interatomic_distances[atom_key])[np.array(store_donor_atoms)])
                    # print(store_donor_atoms)
                    dAtoms_dict = {}
                    for stored_donor_atom in store_donor_atoms:
                        dAtoms_dict[stored_donor_atom] = np.array(interatomic_distances[atom_key])[stored_donor_atom]
                    # print(dAtoms_dict)
                    sorted_dAtoms_dict = {k: v for k, v in sorted(dAtoms_dict.items(), key=lambda item: item[1])}
                    bidentate_indices.extend([list(sorted_dAtoms_dict.keys())[0], list(sorted_dAtoms_dict.keys())[1]])
                    ligand_start_idx = bidentate_indices[1:3]
                    # print(bidentate_indices)
            break

    store_atoms = []
    store_atoms.extend(ligand_start_idx)
    mol_graph = {}

    for (atom_ligand, atom_type) in zip(ligand_start_idx, elements[ligand_start_idx]):      
        if (atom_type in list(atom_covalent_max_dist.keys())) and (atom_type in donor_atoms):
            nr_of_bonds = atom_covalent_max_dist[atom_type] 
            # print(atom_type)
            idx = np.where(np.array(interatomic_distances[atom_ligand]) <= atom_covalent_max_dist[atom_type])
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
            nr_of_bonds = atom_covalent_max_dist[atom_type] 
             
            idx = np.where(np.array(interatomic_distances[atom_key]) <= atom_covalent_max_dist[atom_type])
            idx = idx[:len(idx) + 1]

            # bond_distances = np.array(interatomic_distances[atom_ligand])[idx]                           
            bonded_atoms = list(np.setdiff1d(idx, np.array([atom_key])))
            
            # for id_atom, atom in enumerate(bonded_atoms):
            #     if atom in store_atoms:
            #         del bonded_atoms[id_atom]
                    
            store_atoms.extend(bonded_atoms)
            mol_graph[atom_key] = bonded_atoms
            
    mol_graph[metal_key] = ligand_start_idx

    # ### Adjacency Matrix (fingerprint)
    # keys=sorted(mol_graph.keys())
    # size=len(keys)

    # adj_matrix = [ [0]*size for i in range(size) ]

    # for a,b in [(keys.index(a), keys.index(b)) for a, row in mol_graph.items() for b in row]:
    #     adj_matrix[a][b] = 1
    # plt.spy(adj_matrix)
    # plt.savefig('spy.png')  
    
    ### Breadth first search algorithm (go through all neighbours and store)
    ligands_atoms_idx = {}
    checklist = []
    try:
        
        if geom == 'OH':
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
                    if ligand_index1 != ligand_index2:
                        if ligand_index1 in list(ligands_atoms_idx[ligand_index2]):               
                            if ligand_index1 not in bidentate_indices:
                                bidentate_indices.append(ligand_index1)
                            if ligand_index2 not in bidentate_indices:
                                bidentate_indices.append(ligand_index2)
                                break
        
    # for elem in ligands_atoms_idx:
    #     print(elements[elem]) 
    except Exception:            
        # ligands_atoms_idx = None      
        print("failed")
    return ligands_atoms_idx,bidentate_indices


def BFS_SP(graph, start, goal):
    explored = []
     
    # Queue for traversing the
    # graph in the BFS
    queue = [[start]]
     
    # If the desired node is
    # reached
    if start == goal:
        print("Same Node")
        return
     
    # Loop to traverse the graph
    # with the help of the queue
    while queue:
        path = queue.pop(0)
        node = path[-1]
         
        # Condition to check if the
        # current node is not visited
        if node not in explored:
            neighbours = graph[node]
             
            # Loop to iterate over the
            # neighbours of the node
            for neighbour in neighbours:
                new_path = list(path)
                new_path.append(neighbour)
                queue.append(new_path)
                 
                # Condition to check if the
                # neighbour node is the goal
                if neighbour == goal:
                    print("Shortest path = ", *new_path)
                    return
            explored.append(node)
 
    # Condition when the nodes
    # are not connected
    print("So sorry, but a connecting"\
                "path doesn't exist :(")
    return new_path


def P_P_bridge_path():
    elem, coord = read_xyz('xtbopt.xyz')
    
    _ , bi_PP ,mol_graph = molecular_graph(elem, coord)    
    
    P_P_bridge = BFS_SP(mol_graph, bi_PP[1], bi_PP[2])
    
    return P_P_bridge
    

# # ligands_indices, bidentate = molecular_graph('xtbopt.xyz')
# # print(bidentate, ligands_indices)
# P_P_bridge_path()
