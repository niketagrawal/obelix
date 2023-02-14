# import numpy as np
# import matplotlib.pyplot as plt

# x = np.linspace(0, 1, 200)
# y = np.linspace(1, 0, 200)

# xArray, yArray = np.meshgrid(x, y)
# plotArray = np.sqrt(xArray**2 + yArray**2)

# fig = plt.figure()
# ax = fig.add_subplot(111)

# ax.set_xticks(range(1), ["S1"])

# ax.imshow(plotArray,
#           cmap=plt.cm.Greens,
#           vmin=0,
#           vmax=1)

# plt.savefig("gradient.png")

from run_workflow import *


def remove_ferrocene_substrate(path_to_ferrocene, substrate_position):

    ferrocenes = glob.glob(os.path.join(path_to_ferrocene,"*"))
    dictionary = {}
    
    for ferrocene in ferrocenes:
        print(ferrocene)
        with open(ferrocene + "/xtbopt.xyz", "r+") as f:
            lines = f.readlines()
            lines[0] = str(int(lines[0]) - len(substrate_position) + 1) + "\n"
            del lines[substrate_position[0]:substrate_position[-1]]

            f.seek(0)
            f.truncate()
            f.writelines(lines)
            f.close()
            elements, coordinates = read_xyz(ferrocene + "/xtbopt.xyz")
            dictionary[os.path.basename(ferrocene)] = ConeAngle(elements, coordinates, find_bidentate(ferrocene + "/xtbopt.xyz")[1]).cone_angle
    
    df = dataframe_from_dictionary(dictionary)
    df.to_excel("ferrocene_cone_angle1.xlsx")
    

path = os.path.join(os.getcwd(), "Workflow", "ferrocenes")    
remove_ferrocene_substrate(path, list(range(30, 40)))