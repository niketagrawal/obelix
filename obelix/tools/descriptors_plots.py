#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import scipy
# from umap import umap 

plt.close('all')
bd = pd.read_excel('descriptors_SP.xlsx').dropna()
oh = pd.read_excel('descriptors_OH.xlsx').dropna()
mace_averaged_bd =bd.groupby(['Cas']).mean()
mace_averaged_oh =oh.groupby(['Cas']).mean()


# plt.figure()
# plt.scatter(mace_averaged_bd['bite_angle'], mace_averaged_oh['bite_angle'])
# a = np.polyfit(mace_averaged_bd['bite_angle'], mace_averaged_oh['bite_angle'], deg = 1)

# a1 = a[0]*np.min(mace_averaged_bd['bite_angle']) + a[1]
# a2 = a[0]*np.max(mace_averaged_bd['bite_angle']) + a[1]
# plt.plot([np.min(mace_averaged_bd['bite_angle']), np.max(mace_averaged_bd['bite_angle'])], [a1, a2], '-')
# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(mace_averaged_bd['bite_angle'], mace_averaged_oh['bite_angle'])
# plt.legend([r'$R^2 = $' + str(np.round(r_value, 2))])
# plt.xlabel('Bite angle BD')
# plt.ylabel('Bite angle OH')


# mace_averaged_bd =  mace_averaged_bd.add_suffix('_bd')    
# print(mace_averaged_bd['bite_angle_bd'])
# # print(mace_averaged_oh['bite_angle'])

# # mace_averaged_oh = mace_averaged_oh.loc[mace_averaged_oh['yield'] > 0.6 ]
# # mace_averaged_oh = mace_averaged_oh.loc[mace_averaged_oh['ee'] > 0 ]
# plt.figure()
# print(len(mace_averaged_oh))
# ee_y = np.array([np.abs(mace_averaged_oh['ee']),mace_averaged_oh['yield']]).T
# k_means = KMeans(n_clusters=2)
# label = k_means.fit_predict(ee_y)
# u_labels = np.unique(label)
# mace_averaged_oh['label'] = label
# plt.scatter(ee_y[:,1], ee_y[:,0], c=label, s=25)
# plt.savefig('ee_y.png',  dpi = 400)

# mace_averaged_oh.to_excel('new_oh_kmeans.xlsx')

#%%
import matplotlib.pyplot as plt
plt.close('all')
import pandas as pd
import numpy as np
from pandas.plotting import table 
from umap import UMAP
from sklearn.decomposition import PCA

df = mace_averaged_oh
geom_keys = ['bite_angle', 'cone_angle']
steric_keys = ['buried_volume', 'dispersion', 'sasa', 'dipole']
elec_keys = ['HOMO_LUMO_gap', 'ea', 'ip', 'nucleofugality', 'electrofugality', 'electrophilicity', 'nucleophilicity']
# df = 0-df.loc()


pca = PCA(n_components=1)
u_geom = pca.fit_transform(df[geom_keys]) 
u_steric = pca.fit_transform(df[steric_keys])
u_elec = pca.fit_transform(df[elec_keys])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Gen colors 
# reds = np.zeros((len(df['yield']), 3))
# for i in range(len(df['yield'])):
#     reds[i,:] = [df['yield'][i],  df['yield'][i]/10, df['yield'][i]/10]

plt.rcParams['image.cmap'] = 'RdBu'
print(max(np.abs(df['ee'])))
p = ax.scatter(u_geom, u_elec, u_steric, c=np.abs(df['yield']), s=400)
ax.set_xlabel('Geometric') # geom effects
ax.set_ylabel('Electronic') 
ax.set_zlabel('Steric')
fig.colorbar(p, fraction=0.025, pad=0.05)


fig = plt.figure()
plt.xlabel('Geometric')
plt.ylabel('Electronic')
p = plt.scatter(u_geom, u_elec, c=np.abs(df['yield']), s=100)
fig.colorbar(p)

fig = plt.figure()
plt.xlabel('Geometric')
plt.ylabel('Steric')
p = plt.scatter(u_geom, u_elec, c=np.abs(df['yield']), s=100)                     
fig.colorbar(p)

fig = plt.figure()
plt.xlabel('Electronic')
plt.ylabel('Steric')
p = plt.scatter(u_elec, u_steric, c=np.abs(df['yield']), s=100)
fig.colorbar(p)

# plt.savefig('Interpretable_MAP_PCA.png', dpi = 400, bbox_inches='tight')
fig = plt.figure()
plt.xlabel('Electronic')
plt.ylabel('Yield')
p = plt.scatter(u_elec, df['yield'], c=np.abs(df['label']), s=100)
fig.colorbar(p)

fig = plt.figure()
plt.xlabel('Steric')
plt.ylabel('Yield')
p = plt.scatter(u_steric, df['yield'], c=np.abs(df['label']), s=100)
fig.colorbar(p)

df['Cas'] = df.index
dataframe = pd.DataFrame()
dataframe['steric'] = list(u_steric.reshape(45, -1))
dataframe['elec'] = list(u_elec.reshape(45,-1))
dataframe['geom'] = list(u_geom.reshape(45, -1))
dataframe['yield'] = df['yield'].to_list()
dataframe['Cas'] = df['Cas'].to_list()
dataframe['ee'] = df['ee'].to_list()

dataframe.to_excel('PCA2.xlsx')