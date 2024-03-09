# microglia clusters 0&1
# library(Seurat)
# m = readRDS('/home/brasel/SingleCellProjects/dataObjects/microglia.rds')
# m = subset(m,subset=seurat_clusters %in% c(0,1))
# m = GetAssayData(m,slot='counts')
# write.csv(m,file = '/home/brasel/SingleCellProjects/dataObjects/mic_01.csv')

from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from tqdm import tqdm
import pickle as pkl
import pandas as pd
import numpy as np
import os

from velocyto.estimation import colDeltaCorpartial
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from scipy.stats import norm as normal

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib
import logging
import sys
from typing import List, Optional

RF_KWARGS = {
    'n_jobs': 1,
    'n_estimators': 1000,
    'max_features': 'sqrt'
}

GBM_KWARGS = {
    'learning_rate': 0.01,
    'n_estimators': 500,
    'max_features': 0.1
}

SKLEARN_REGRESSOR_FACTORY = {
    'RF': RandomForestRegressor,
    'GBM': GradientBoostingRegressor
}

DEFAULT_REGRESSOR_PARAMS = {
    'RF': RF_KWARGS,
    'GBM': GBM_KWARGS
}

with open('snakemake_TF_plot.pkl','wb') as f:
    pkl.dump(snakemake,f)

def calculate_UMAP_embeddings(
        obj = None,
        perturbed_matrix = None,
        label = None):
    o = np.log(obj+1)
    pm = np.log(perturbed_matrix+1)
    #R
    #m = readRDS('/home/brasel/SingleCellProjects/dataObjects/microglia.rds')
    #redeuc = m[['umap']]
    #data.use = GetAssayData(m,slot='data')
    #cell.embeddings <- Embeddings(object = redeuc)
    #new.feature.loadings.full <- data.use %*% cell.embeddings
    umap_loadings = pd.read_csv(str(snakemake.input['feature_loadings']),index_col=0)
    umap_loadings = umap_loadings.loc[[g for g in umap_loadings.index.values if g in obj.columns.values],:]
    umap = o[umap_loadings.index.values].dot(umap_loadings)
    return(umap)


def permute_rows_nsign(A: np.ndarray) -> None:
    """Permute in place the entries and randomly switch the sign for each row of a matrix independently.
    From celloracle
    """
    plmi = np.array([+1, -1])
    for i in range(A.shape[0]):
        np.random.shuffle(A[i, :])
        A[i, :] = A[i, :] * np.random.choice(plmi, size=A.shape[1])

def _project_perturbation_in_embedding(
    original_matrix, 
    perturbed_matrix, 
    embedding, 
    save,
    sigma_corr = 0.05, n_cpu = 1):
    #based on celloracle/velocyto code
    #delta_matrix = perturbed_matrix.to_numpy().astype('double') - original_matrix.to_numpy().astype('double')
    #embedding = embedding.to_numpy()
    delta_matrix = perturbed_matrix.to_numpy().astype('double') - original_matrix.to_numpy().astype('double')
    delta_matrix_random =  delta_matrix.copy()
    permute_rows_nsign(delta_matrix_random)

    n_neighbors = int(perturbed_matrix.shape[0] / 5) #default from cell oracle
    nn = NearestNeighbors(n_neighbors = n_neighbors + 1, n_jobs = n_cpu)
    nn.fit(embedding)
    embedding_knn = nn.kneighbors_graph(mode = 'connectivity')

    # Pick random neighbours and prune the rest
    neigh_ixs = embedding_knn.indices.reshape((-1, n_neighbors + 1))
    p = np.linspace(0.5, 0.1, neigh_ixs.shape[1])
    p = p / p.sum()

    # There was a problem of API consistency because the random.choice can pick the diagonal value (or not)
    # resulting self.corrcoeff with different number of nonzero entry per row.
    # Not updated yet not to break previous analyses
    # Fix is substituting below `neigh_ixs.shape[1]` with `np.arange(1,neigh_ixs.shape[1]-1)`
    # I change it here since I am doing some breaking changes
    sampling_ixs = np.stack([np.random.choice(neigh_ixs.shape[1],
                                            size=(int(0.3 * (n_neighbors + 1)),),
                                            replace=False,
                                            p=p) for i in range(neigh_ixs.shape[0])], 0)
    neigh_ixs = neigh_ixs[np.arange(neigh_ixs.shape[0])[:, None], sampling_ixs]
    nonzero = neigh_ixs.shape[0] * neigh_ixs.shape[1]
    embedding_knn = sparse.csr_matrix((np.ones(nonzero),
                                        neigh_ixs.ravel(),
                                        np.arange(0, nonzero + 1, neigh_ixs.shape[1])),
                                        shape=(neigh_ixs.shape[0],
                                                neigh_ixs.shape[0]))

    corrcoef = colDeltaCorpartial(perturbed_matrix.T.astype('double'), delta_matrix.T.astype('double'), neigh_ixs,  threads = n_cpu)
    corrcoef[np.isnan(corrcoef)] = 1
    np.fill_diagonal(corrcoef, 0)

    transition_prob = np.exp(corrcoef / sigma_corr) * embedding_knn.A
    transition_prob /= transition_prob.sum(1)[:, None]
    
    unitary_vectors = embedding.T[:, None, :] - embedding.T[:, :, None]  # shape (2,ncells,ncells)
    unitary_vectors /= np.linalg.norm(unitary_vectors, ord=2, axis=0)  # divide by L2
    np.fill_diagonal(unitary_vectors[0, ...], 0)  # fix nans
    np.fill_diagonal(unitary_vectors[1, ...], 0)

    delta_embedding = (transition_prob * unitary_vectors).sum(2)
    delta_embedding = delta_embedding - ((embedding_knn.A * unitary_vectors).sum(2) / embedding_knn.sum(1).A.T)
    delta_embedding = delta_embedding.T
    with open(save, "wb") as f:
        pkl.dump(delta_embedding, f)
    return delta_embedding

def _calculate_grid_arrows(embedding, delta_embedding, offset_frac, n_grid_cols, n_grid_rows, n_neighbors, n_cpu):
    #prepare grid
    min_x = min(embedding[:, 0])
    max_x = max(embedding[:, 0])
    min_y = min(embedding[:, 1])
    max_y = max(embedding[:, 1])
    offset_x = (max_x - min_x) * offset_frac
    offset_y = (max_y - min_y) * offset_frac
    #calculate number of points underneath grid points
    x_dist_between_points = (max_x - min_x) / n_grid_cols
    y_dist_between_points = (max_y - min_y) / n_grid_rows
    minimal_distance = np.mean([y_dist_between_points, x_dist_between_points]) #will be used to mask certain points in the grid

    grid_x, grid_y = np.meshgrid(
        np.linspace(min_x + offset_x, max_x - offset_x, n_grid_cols),
        np.linspace(min_y + offset_y, max_y - offset_y, n_grid_rows)
    )
    grid_xy = np.array([np.hstack(grid_x), np.hstack(grid_y)]).T

    #find neighbors of gridpoints
    nn = NearestNeighbors(n_neighbors = n_neighbors, n_jobs = n_cpu)
    nn.fit(embedding)
    dists, neighs = nn.kneighbors(grid_xy)

    std = np.mean([abs(g[1] - g[0]) for g in grid_xy])
    # isotropic gaussian kernel
    gaussian_w = normal.pdf(loc=0, scale=0.5*std, x=dists)
    total_p_mass = gaussian_w.sum(1)

    uv = (delta_embedding[neighs] * gaussian_w[:, :, None]).sum(1) / np.maximum(1, total_p_mass)[:, None]

    #mask points in the grid which don't have points of the embedding underneath them
    mask = dists.min(1) < minimal_distance

    return grid_xy, uv, mask

def plot_perturbation_effect_in_embedding(
    delta_embedding,
    embedding,
    fullEmbedding,
    clusters,
    grid_offset_frac: Optional[float] = 0.005,
    grid_n_cols: Optional[int] = 25,
    grid_n_rows: Optional[int] = 25,
    grid_n_neighbors: Optional[int] = 25,
    n_cpu: Optional[int] = 1,
    figsize: Optional[tuple] = (6.4, 6.4),
    title = None,
    save: Optional[str] = None):
    #embedding = embedding.to_numpy()
    grid_xy, uv, mask = _calculate_grid_arrows(
        embedding=embedding, 
        delta_embedding=delta_embedding,
        offset_frac=grid_offset_frac,
        n_grid_cols=grid_n_cols,
        n_grid_rows=grid_n_rows,
        n_neighbors=grid_n_neighbors,
        n_cpu=n_cpu)
    distances = np.sqrt((uv**2).sum(1))
    norm = matplotlib.colors.Normalize(vmin=0.15, vmax=0.5, clip=True)
    scale = lambda X: [(x - min(X)) / (max(X) - min(X)) for x in X]
    uv[np.logical_or(~mask, np.array(scale(distances)) < 0.15)] = np.nan
    #log.info('Plotting')

    # Replace streamplot with quiver
    X = grid_xy.reshape(grid_n_cols,grid_n_rows, 2)[:, :, 0]
    Y = grid_xy.reshape(grid_n_cols,grid_n_rows, 2)[:, :, 1]
    U = uv.reshape(grid_n_cols,grid_n_rows, 2)[:, :, 0]
    V = uv.reshape(grid_n_cols,grid_n_rows, 2)[:, :, 1]
    C = np.array(scale(distances)).reshape(grid_n_cols, grid_n_rows)  # color based on distances

    # Define colors
    if snakemake.params['cellType'] == 'mic':
        colors = ["#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3"]
    else:
        colors = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]
    cmap = mcolors.ListedColormap(colors)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(fullEmbedding[:, 0], fullEmbedding[:, 1], s=10, c=cluster['x'], cmap=cmap, alpha=0.4, edgecolor=(0,0,0,1), lw= 0.3)

    ax.quiver(
        X, Y, U, V,
        C,
        cmap='Greys',
        zorder=10,
        norm=norm,
        width=0.005  # adjust as needed
    )
    #ax.streamplot(
    #        grid_xy.reshape(grid_n_cols,grid_n_rows, 2)[:, :, 0],
    #        grid_xy.reshape(grid_n_cols,grid_n_rows, 2)[:, :, 1],
    #        uv.reshape(grid_n_cols,grid_n_rows, 2)[:, :, 0],
    #        uv.reshape(grid_n_cols,grid_n_rows, 2)[:, :, 1], 
    #        density = 3, 
    #        color = 'red',#np.array(scale(distances)).reshape(grid_n_cols, grid_n_rows),
    #        cmap = 'Greys', 
    #        zorder = 10, 
    #        norm = norm,
    #        linewidth = 0.5)
    ax.set_title(title)
    if save is not None:
        fig.savefig(save)
    else:
        plt.show(fig)
        return ax

#map perturbed expression to original microglia PCA/UMAP space
# https://www.notion.so/loganlabbook/Sex-Differences-d379f0f2cd624cc0a505c12895158934?pvs=4#88f6d492b4f34af5a64a1fbf34b07317
# https://github.com/satijalab/seurat/blob/ca4c48b6300d7c10857e2d59a8ddee2858c7e2fc/R/dimensional_reduction.R#L282
# m = ProjectDim(m,reduction = 'umap')

print('loading data...')
obj = pd.read_csv(snakemake.input['cellStateMat'],index_col=0)
obj = obj.T

#fullEmbedding = pd.read_csv('microglia/data/umap_embeddings_all_mic.csv',index_col=0)
fullEmbedding = np.genfromtxt(str(snakemake.input['fullEmbed']), delimiter=',', skip_header=1, usecols=(1,2))
pertEmbedding= np.genfromtxt(str(snakemake.input['pertEmbed']), delimiter=',', skip_header=1, usecols=(1,2))
cluster = pd.read_csv(str(snakemake.input['clusters']), index_col=0)
print('data loaded')

with open(snakemake.input['perturbed_matrix'], "rb") as f:
    perturbed_matrix = pkl.load(f)

delta_file = snakemake.output['delta_matrix']
print('calculating delta embedding')
if os.path.isfile(delta_file):
    with open(delta_file, "rb") as f:
        delta_embedding = pkl.load(f)
else:
    delta_embedding = _project_perturbation_in_embedding(np.log(obj+1), perturbed_matrix, pertEmbedding, save = delta_file)

print('plotting')
title = f"{snakemake.params['cellType']}.{snakemake.params['cellState']}_{snakemake.params['TF']}"
plot_perturbation_effect_in_embedding(delta_embedding, pertEmbedding, fullEmbedding, clusters = cluster, title = title, save = snakemake.output['plot'])
