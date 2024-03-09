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
from scipy.stats import norm as normal

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

with open ("test.pkl", "wb") as f:
    pkl.dump(snakemake, f)

def _do_one_round_of_simulation(original_matrix, perturbed_matrix, regressors, smpl_counts):
    new_exp_mtx = perturbed_matrix.copy()
    original_matrix['logCounts'] = np.log([smpl_counts[s.split('_')[0]] for s in original_matrix.index.values]) ########################
    perturbed_matrix['logCounts'] = np.log([smpl_counts[s.split('_')[0]] for s in perturbed_matrix.index.values]) ########################
    for gene in tqdm(regressors.keys(), total = len(regressors.keys()), leave = False):
        gene_regressor = regressors[gene][-1]
        TF_original_exp = original_matrix[regressors[gene][0: len(regressors[gene]) -1]+['logCounts']].to_numpy() ########################
        TF_perturbed_exp = perturbed_matrix[regressors[gene][0: len(regressors[gene]) -1]+['logCounts']].to_numpy() ########################
        if not all((TF_original_exp == TF_perturbed_exp).ravel()):
            gene_predicted_exp_orig = gene_regressor.predict(TF_original_exp)
            gene_predicted_exp_perturbed = gene_regressor.predict(TF_perturbed_exp)
            fc = gene_predicted_exp_perturbed / gene_predicted_exp_orig
            new_exp_mtx[gene] = new_exp_mtx[gene] * fc
        else:
            continue
    return new_exp_mtx

def getSampleCounts(obj):
    smpl = [s.split('_')[0] for s in obj.index.values]
    smpl_counts = {}
    for s in smpl:
        if s not in smpl_counts.keys():
            smpl_counts[s] = 0
        smpl_counts[s] += 1
    return(smpl_counts)

#simulate perturbation
def simulate_perturbation(
    obj = None,
    perturbation: dict = {'CEBPB': 0},
    regressors: dict = None,
    keep_intermediate = False,
    save = None,
    save_mtx = None,
    n_iter = 5):
    if save is None:
        raise ValueError(f'save must be provided to create unique file names for each model')
    #initialize
    if keep_intermediate:
            perturbation_over_iter = {}
    smpl_counts = getSampleCounts(obj) ########################
    original_matrix = np.log(obj.copy() + 1) ########################
    perturbed_matrix = original_matrix.copy()
    for gene in perturbation.keys():
        perturbed_matrix[gene] = perturbed_matrix[gene] * perturbation[gene]
    if keep_intermediate:
        perturbation_over_iter['0'] = original_matrix
        perturbation_over_iter['1'] = perturbed_matrix
    #do several iterations of perturbation
    for i in range(n_iter):
        new_matrix = _do_one_round_of_simulation(original_matrix, perturbed_matrix, regressors, smpl_counts)
        original_matrix = perturbed_matrix.copy()
        perturbed_matrix = new_matrix.copy()
        if keep_intermediate:
            perturbation_over_iter[str(i + 2)] = perturbed_matrix
    #save simulated perturbation pickle
    with open(save, "wb") as f:
        pkl.dump(perturbed_matrix, f)
    return perturbed_matrix

#map perturbed expression to original microglia PCA/UMAP space
# https://www.notion.so/loganlabbook/Sex-Differences-d379f0f2cd624cc0a505c12895158934?pvs=4#88f6d492b4f34af5a64a1fbf34b07317
# https://github.com/satijalab/seurat/blob/ca4c48b6300d7c10857e2d59a8ddee2858c7e2fc/R/dimensional_reduction.R#L282
# m = ProjectDim(m,reduction = 'umap')

print('loading data...')
obj = pd.read_csv(snakemake.input['cellStateMat'],index_col=0)
obj = obj.T

with open(snakemake.input['regressors'], "rb") as f:
    regressors = pkl.load(f)
print('data loaded')

pm_file = snakemake.output['perturbed_matrix']
print('simulating perturbation')
if os.path.isfile(pm_file):
    with open(pm_file, "rb") as f:
        perturbed_matrix = pkl.load(f)
else:
    perturbed_matrix = simulate_perturbation( obj = obj, perturbation = {snakemake.params['TF']: 0}, regressors = regressors, save = pm_file)
#perturbed_matrix = simulate_perturbation( obj = obj, perturbation = {snakemake.params['TF']: 0}, regressors = regressors, save = snakemake.output['perturbed_matrix'])
