#%%
import numpy as np
import morfeus as mf
import os

def dihedral(p):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

Ir = np.array([-1.25785409922576, -1.20776561791263, 1.20123689358246]) # N
P1 = np.array([-2.25983878035412, 0.63938236826987, 0.02615142859852]) # CA
P2 = np.array([0.83234347300909, -0.99047236842121, 0.16680357180095]) # C
N = np.array([-0.51778338173685, -0.06749071682697, 2.56048254471549]) # O

elements, coordinates = mf.read_xyz(os.getcwd() + '/new_descriptors/xtbopt.xyz')

elem_coords = []
for elem_index, element in enumerate(elements):
    if element in ['P', 'N', 'Ir']:
        elem_coords.append(coordinates[elem_index])
    
dihedral(elem_coords)

# print(elements, coordinates)
# print('Dihedral_test', 'dihedral_P_P_sub_C', dihedral(np.array([] )))

