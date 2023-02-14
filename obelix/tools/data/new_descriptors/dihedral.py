#%%
import numpy as np
import morfeus as mf
import os, glob
import pandas as pd

def dataframe_from_dictionary(dictionary):
    Dataframe = pd.DataFrame.from_dict({i: dictionary[i]
                           for i in dictionary.keys()}, orient='index')
    return Dataframe


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


ferro = glob.glob(os.getcwd() + "/ferro/*")
# elements, coordinates = mf.read_xyz(os.getcwd() + '/xtbopt.xyz')


# 13-14,20,22,27-28


dihedral_indices = [[22, 14, 13, 20],
                    [27, 20, 13, 14],
                    [28, 27, 20, 13],
                    [28, 22, 14, 13]]

dihedral_ = {}

for ferrocene in ferro:

    elements, coordinates = mf.read_xyz(ferrocene + "/xtbopt.xyz")
    # for element in np.array([13, 14, 20, 22, 27, 28]) - 1:
    #     elem_coords.append(coordinates[element])
    
    elem_coords_all = []
    i = 0
    
    dihedral_angles = {}
    for dihedral_angle in dihedral_indices:
        elem_coords = []
        for atom in dihedral_angle:
            elem_coords.append(coordinates[atom])
        # elem_coords_all.append(elem_coords)
        i+=1 
        dihedral_angles[i] = dihedral(elem_coords)
    dihedral_[os.path.basename(ferrocene)] = dihedral_angles    
        
df = dataframe_from_dictionary(dihedral_)
df.to_excel("dihedral.xlsx")

# Calculate 4 dihedral-angles

#1 22 - 14 - 13 - 20 P  - (CFe1) - (CFe2) - (C-Pc)
#2 27 - 20 - 13 - 14 Pc - (C-Pc) - (CFe1) - (CFe2)
#3 28 - 27 - 20 - 13 Ir - (Pc) - (C-Pc) - (CFe1)
#4 28 - 22 - 14 - 13 Ir - P - (CFe2) - (CFe1)

# 27-28 - p-c, Ir
# 22 - p
# 13-14 -> Fe-C, Fe-C-(P)
# 20 - Carbon bridge to ferrocene


# %%
