import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from scipy.io import mmread
import pickle as pkl
import matplotlib.pyplot as plt

# Load the GEM and UMAP data from file
geneFilt = pd.read_csv(str(snakemake.input['geneFilt']),header=None)
umat = pd.read_csv(str(snakemake.input['umap']), index_col=0)
sparse_matrix = mmread(str(snakemake.input['cellMat']))
allGenes = pd.read_csv(str(snakemake.input['allGenes']))
gem = pd.DataFrame(np.log(sparse_matrix.T.toarray()+1))
#del sparse_matrix
gem.index = umat.index
gem.columns = allGenes['x']
gem = gem.loc[:,geneFilt[0].values]

# Split the GEM and UMAP data into train, test, and validation sets based on the cluster vector
cluster = pd.read_csv(str(snakemake.input['clusters']), index_col=0)
clusterOneHot = pd.get_dummies(cluster, columns=['x']).astype(float)
gem = pd.concat([gem, clusterOneHot], axis=1)

# Load Scaler
with open(snakemake.input['scaler'], 'rb') as f:
    scaler = pkl.load(f)

# Load model
model = tf.keras.models.load_model(str(snakemake.input['model']))

# Compare the predicted UMAT values to the actual UMAT values by plotting them
full_umat = model.predict(scaler.transform(gem))

#subset data
umat_subset = full_umat[np.where(cluster['x'] == snakemake.params['cellState'])[0]]

#plot the underlying full UMAP
plt.scatter(full_umat[:,0], full_umat[:,1], c=cluster['x'],s = 1)
#plot the original cell state on top
plt.scatter(umat_subset[:,0], umat_subset[:,1], c='blue',s = 1)

# Use the trained model to predict UMAT values for new GEM data
with open(snakemake.input['pertMat'], 'rb') as f:
    pert_gem = pkl.load(f) # already log transformed
#make cluster one hot. set to cluster 0
pert_gem = pd.concat([pert_gem, clusterOneHot.loc[pert_gem.index.values]], axis=1)
new_umat = model.predict(scaler.transform(np.log(pert_gem+1)))

#plot the perturbation on top
plt.scatter(new_umat[:,0], new_umat[:,1], c='red',s = 1)
plt.savefig(snakemake.output['plot'])
plt.clf()

# Save new_umat to csv file
new_umat = pd.DataFrame(new_umat, index=pert_gem.index)
new_umat.to_csv(snakemake.output['perturbed_cell_embeddings'])