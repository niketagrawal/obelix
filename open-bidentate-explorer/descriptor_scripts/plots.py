# Find descriptor corr.
# Build dataframe where [SP, OH]

import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
# from sklearn import linear_model
# from sklearn.metrics import r2_score
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


#####################################################################################

df_SP =pd.DataFrame(pd.read_excel("Build_for_py.xlsx", "SP_COMMON"))
df_OH =pd.DataFrame(pd.read_excel("Build_for_py.xlsx", "OH_COMMON"))

boltzmann_avg_cols = [col for col in df_SP if 'boltzmann' in col]
nr_of_desc = len(boltzmann_avg_cols)
combs_2_descriptors = int(math.factorial(nr_of_desc)/math.factorial(nr_of_desc - 2))

nr_of_plot_col = 5
nr_of_plot_row = (combs_2_descriptors + nr_of_desc)/nr_of_plot_col
fig, ax = plt.subplots(int(nr_of_plot_row), nr_of_plot_col)

trace_duplicates = []


for i, _ in enumerate(boltzmann_avg_cols):
    for j, _ in enumerate(boltzmann_avg_cols):
        trace_duplicates.append([i, j])

row = 0
col = 0

# for [i, j] in trace_duplicates:    
    
#     ### LINEAR REGRESSION    

#     X = df_SP[boltzmann_avg_cols[i]].values
#     y = df_OH[boltzmann_avg_cols[j]].values
    
#     X = np.vstack([X, np.ones(len(X))]).T
#     x_min = np.min(X)
#     x_max = np.max(X)
#     model, resid = np.linalg.lstsq(X, y, rcond=None)[:2]
    
#     R2 = 1 - resid/y.size/y.var()
#     # mx + c 
#     m = model[0]
#     c = model[1]
    
#     # R2 = r2_score(df_SP[boltzmann_avg_cols[i]], predict_data)
#     ### SAVE SEP
    
#     plt.figure()
#     plt.scatter(df_SP[boltzmann_avg_cols[i]], df_OH[boltzmann_avg_cols[j]])
#     plt.xlabel(boltzmann_avg_cols[i])
#     plt.ylabel(boltzmann_avg_cols[j])
#     plt.plot([x_min, x_max], [model[0]*x_min + model[1], model[0]*x_max + model[1]])
#     plt.legend(['{}*x + {}'.format(np.round(m, 2), np.round(c)), r'$R^2 = {}$'.format(R2)]) 
#     # ## SAVE ALL IN ONE
#     # ax[row][col].scatter(df_SP[boltzmann_avg_cols[i]], df_OH[boltzmann_avg_cols[j]])
#     # ax[row][col].set_xlabel(boltzmann_avg_cols[i])
#     # ax[row][col].set_ylabel(boltzmann_avg_cols[j])
#     # col+=1
    
#     # if col == 5:
#     #     row+=1
#     #     col=0
    
#     plt.savefig('figures/corr_SP_{}_OH_{}.png'.format(boltzmann_avg_cols[i], boltzmann_avg_cols[j]),figsize=(32, 10), bbox_inches='tight')


# Clustering of linear models.


# for [i, j] in trace_duplicates:
    
#     # Build dataframe for KMeans
    
#     x = df_SP[boltzmann_avg_cols[i]].values
#     y = df_OH[boltzmann_avg_cols[j]].values
#     X = np.vstack((x, y)).T
#     range_n_clusters = [2, 3, 4]

#     for n_clusters in range_n_clusters:
#         # Create a subplot with 1 row and 2 columns
#         fig, (ax1, ax2) = plt.subplots(1, 2)
#         fig.set_size_inches(18, 7)

#         # The 1st subplot is the silhouette plot
#         # The silhouette coefficient can range from -1, 1 but in this example all
#         # lie within [-0.1, 1]
#         ax1.set_xlim([-0.1, 1])
#         # The (n_clusters+1)*10 is for inserting blank space between silhouette
#         # plots of individual clusters, to demarcate them clearly.
#         ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

#         # Initialize the clusterer with n_clusters value and a random generator
#         # seed of 10 for reproducibility.
#         clusterer = KMeans(n_clusters=n_clusters, random_state=10)
#         cluster_labels = clusterer.fit_predict(X)

#         # The silhouette_score gives the average value for all the samples.
#         # This gives a perspective into the density and separation of the formed
#         # clusters
#         silhouette_avg = silhouette_score(X, cluster_labels)
#         print(
#             "For n_clusters =",
#             n_clusters,
#             "The average silhouette_score is :",
#             silhouette_avg,
#         )

#         # Compute the silhouette scores for each sample
#         sample_silhouette_values = silhouette_samples(X, cluster_labels)

#         y_lower = 10
#         for n in range(n_clusters):
#             # Aggregate the silhouette scores for samples belonging to
#             # cluster i, and sort them
#             ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == n]

#             ith_cluster_silhouette_values.sort()

#             size_cluster_i = ith_cluster_silhouette_values.shape[0]
#             y_upper = y_lower + size_cluster_i

#             color = cm.nipy_spectral(float(n) / n_clusters)
#             ax1.fill_betweenx(
#                 np.arange(y_lower, y_upper),
#                 0,
#                 ith_cluster_silhouette_values,
#                 facecolor=color,
#                 edgecolor=color,
#                 alpha=0.7,
#             )

#             # Label the silhouette plots with their cluster numbers at the middle
#             ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(n))

#             # Compute the new y_lower for next plot
#             y_lower = y_upper + 10  # 10 for the 0 samples

#         ax1.set_title("The silhouette plot for the various clusters.")
#         ax1.set_xlabel("The silhouette coefficient values")
#         ax1.set_ylabel("Cluster label")

#         # The vertical line for average silhouette score of all the values
#         ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

#         ax1.set_yticks([])  # Clear the yaxis labels / ticks
#         ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

#         # 2nd Plot showing the actual clusters formed
#         colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
#         ax2.scatter(X[:, 0], X[:, 1], marker=".", s=40, lw=0, alpha=0.7, c=colors, edgecolor="k")

#         # Labeling the clusters
#         centers = clusterer.cluster_centers_

#         # Draw white circles at cluster centers

#         ax2.scatter(
#             centers[:, 0],
#             centers[:, 1],
#             marker="o",
#             c="white",
#             alpha=1,
#             s=200,
#             edgecolor="k",
#         )

#         for l, c in enumerate(centers):
#             ax2.scatter(c[0], c[1], marker="$%d$" % l, alpha=1, s=50, edgecolor="k")

#         ax2.set_title("The visualization of the clustered data.")
#         ax2.set_xlabel(boltzmann_avg_cols[i])
#         ax2.set_ylabel(boltzmann_avg_cols[j])

#         plt.suptitle(
#             "Silhouette analysis for KMeans clustering on sample data with n_clusters = %d"
#             % n_clusters,
#             fontsize=14,
#             fontweight="bold",
#         )

#         plt.savefig('figures/SP_OH_Descriptor_corr_clustering/corr_SP_{}_OH_{}_cluster.png'.format(boltzmann_avg_cols[i], boltzmann_avg_cols[j]),figsize=(32, 10), bbox_inches='tight')


def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')



# df2 = df_SP[df_SP['boltzmann_bite_angle']<70][['sub_1',	'sub_2', 'sub_3' ,'sub_4']]
    
# print_full(x) 

# HOMO_LUMO_GAP_SP = df_SP['boltzmann_homo'] - df_SP['boltzmann_lumo']

# HOMO_LUMO_GAP_OH = df_OH['boltzmann_homo'] - df_OH['boltzmann_lumo']


# plt.figure()
# plt.scatter(np.log(-HOMO_LUMO_GAP_SP), df_OH['boltzmann_energy'])

# plt.savefig('figures/homo_lumo/ln_homo_lumo_SP_energy_OH.png')

dataframe_descriptors_SP = df_SP[boltzmann_avg_cols]
dataframe_descriptors_SP.columns = [str(col) + '_SP' for col in dataframe_descriptors_SP.columns]
dataframe_descriptors_OH = df_OH[boltzmann_avg_cols]
dataframe_descriptors_OH.columns = [str(col) + '_OH' for col in dataframe_descriptors_OH.columns]


concatenated = dataframe_descriptors_SP.join(dataframe_descriptors_OH)
correlation_matrix = concatenated.corr()**2
print(concatenated)
correlation_matrix.to_excel('figures/corr1.xlsx')