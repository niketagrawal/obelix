# from chemplot import Plotter
# import chemplot as cp
from pandas import read_excel
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
from sklearn.cluster import KMeans
from os import mkdir, chdir
from rdkit.Chem import AllChem, MolFromMolFile, MolFromSmiles, MolToMolFile, rdmolfiles
from rdkit.Chem import Draw
from utilities import get_cluster_centroid_coord


def generate_image(cas, mol):
    # mol = MolFromMolFile('x1.mol')
    mol = MolFromMolFile(mol)

    # molFile = MolToMolFile(mol, cas + '.mol')
    
    Draw.MolToFile(mol,"images/{}.png".format(cas),size=(400,400))


generate_image('174467-31-3','174467-31-3.mol')

# # C12C3(P(C4CCCCC4)C4CCCCC4)C4(C5[Fe]1341346(C7C1C3C4C67)C25)P(C1=CC=CC=C1)C1=CC=CC=C1

# ### Difference between covalent SMILES and .struc. SMILES

# # generate_image('Ligand1', 'C12C3(P(C4CCCCC4)C4CCCCC4)C4(C5[Fe]1341346(C7C1C3C4C67)C25)P(C1=CC=CC=C1)C1=CC=CC=C1')


# data = read_excel('ligands.xlsx')

# data_to_array = np.array(data)

# for (cas, smiles) in zip(data['Cas'], data['smiles']):
#     generate_image(cas, smiles)

# # print(my_smiles_list)


# cp =  Plotter.from_smiles(data['smiles'])
# umap_data = cp.umap()
# UMAP1 = np.array(umap_data['UMAP-1'])
# UMAP2 = np.array(umap_data['UMAP-2'])

# data = np.vstack((UMAP1, UMAP2)).T
# k_means = KMeans(n_clusters=16)

# label = k_means.fit_predict(data)
# u_labels = np.unique(label)
# umap_data['label'] = label

# plt.figure()

# plt.scatter(data[:,0], data[:,1], c=k_means.labels_.astype(float), s=1)

# for unique_label in u_labels:
#     plt.scatter(data[label == unique_label, 0], data[label == unique_label, 1], label = 'Cluster {}'.format(unique_label), s=3.8)

# plt.legend(bbox_to_anchor = (-.1, 1))
# plt.savefig('clusters_umap.png', bbox_inches='tight', dpi=760)
# plt.show()
# umap_data.to_excel('umap_clusters.xlsx')



#%% 

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import glob
from PIL import Image
from matplotlib import style
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
style.use('ggplot')
# Generate data x, y for scatter and an array of images.

x = np.array(pd.read_excel('umap_clusters.xlsx')['UMAP-1'])
y = np.array(pd.read_excel('umap_clusters.xlsx')['UMAP-2'])

arr = np.empty((len(x),400,400))
images = glob.glob('images/*.png')
# print(np.shape(arr)[0])
for index in range(np.shape(arr)[0]):
    image = images[index]
    image_array = plt.imread(image)
    image_array = Image.open(images[index]).convert('L')
    arr[index] = image_array

# create figure and plot scatter
fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot(x,y, ls="", marker="o", markersize=1)

ax.set_xlabel('UMAP-1')
ax.set_ylabel('UMAP-2')
ax.set_title('Chemical space clustering')
# create the annotations box
im = OffsetImage(arr[0,:,:], zoom=0.3)
xybox=(50., 50.)
ab = AnnotationBbox(im, (0,0), xybox=xybox, xycoords='data',
        boxcoords="offset points",  pad=0.3,  arrowprops=dict(arrowstyle="->"))
# add it to the axes and make it invisible
ax.add_artist(ab)
ab.set_visible(False)

def hover(event):
    # if the mouse is over the scatter points
    try:
        if line.contains(event)[0]:
            # find out the index within the array from the event

            ind,  = line.contains(event)[1]["ind"]
            # get the figure size
            w,h = fig.get_size_inches()*fig.dpi
            ws = (event.x > w/2.)*-1 + (event.x <= w/2.) 
            hs = (event.y > h/2.)*-1 + (event.y <= h/2.)
            # if event occurs in the top or right quadrant of the figure,
            # change the annotation box position relative to mouse.
            ab.xybox = (xybox[0]*ws, xybox[1]*hs)
            # make annotation box visible
            ab.set_visible(True)
            # place it at the position of the hovered scatter point
            ab.xy =(x[ind], y[ind])
            # set the image corresponding to that point
            im.set_data(arr[ind,:,:])
        else:
            #if the mouse is not over a scatter point8ygfc
            ab.set_visible(False)
        fig.canvas.draw_idle()
    except Exception:
        pass

# add callback for mouse moves
fig.canvas.mpl_connect('motion_notify_event', hover)           
plt.show()