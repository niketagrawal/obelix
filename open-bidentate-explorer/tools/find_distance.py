from morfeus import read_xyz
from numpy.linalg import norm

def find_distance(xyz):
    
    elements, coordinates = read_xyz(xyz)
    
    P_coord = []
    for ind_element, element in enumerate(elements):
        if element == 'P':
            P_coord.append(ind_element)
        if element == 'Ir':
            ind_Ir = ind_element

    d1 = norm(coordinates[ind_Ir]  - coordinates[P_coord[0]])
    d2 = norm(coordinates[ind_Ir]  - coordinates[P_coord[1]])

    return d1, d2
    
    