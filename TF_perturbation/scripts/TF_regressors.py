# microglia clusters 0&1
# library(Seurat)
# m = readRDS('/home/brasel/SingleCellProjects/dataObjects/microglia.rds')
# m = subset(m,subset=seurat_clusters %in% c(0,1))
# m = GetAssayData(m,slot='counts')
# write.csv(m,file = '/home/brasel/SingleCellProjects/dataObjects/mic_01.csv')


from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
import pickle as pkl
import pandas as pd
import numpy as np
import os

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

#with open("/home/brasel/SingleCellProjects/MyProjects/pySCENIC_APOE/microglia/data/TF_hits_mic.pkl", "rb") as f:
with open(str(snakemake.input['TFdict']), "rb") as f:
    genesToTFs = pkl.load(f)
    
def getPredictorTFs(gene):
    #return([g for g in genesToTFs['mic0'][gene].keys() if genesToTFs['mic0'][gene][g] >= 5])
    stateKey = snakemake.params['cellType'] + snakemake.params['cellState']
    return([g for g in genesToTFs[stateKey][gene].keys() if genesToTFs[stateKey][gene][g] >= 5])

def getSampleCounts(obj):
    smpl = [s.split('_')[0] for s in obj.index.values]
    smpl_counts = {}
    for s in smpl:
        if s not in smpl_counts.keys():
            smpl_counts[s] = 0
        smpl_counts[s] += 1
    return(smpl_counts)

#train gene expression models
def train_gene_expression_models(
    obj = None,
    genes = None,#obj.columns.values 
    regressor_type= 'GBM',
    regressor_kwargs= None,
    save = None):
    if regressor_type not in SKLEARN_REGRESSOR_FACTORY.keys():
        raise ValueError(f'Please select a regressor_type from {", ".join(SKLEARN_REGRESSOR_FACTORY.keys())}')
    if regressor_kwargs is None:
        regressor_kwargs = DEFAULT_REGRESSOR_PARAMS[regressor_type]
    #if genes is set to None, predict expression of all genes
    if genes is None:
        genes = obj.columns.values #scplus_obj.gene_names
    smpl_counts = getSampleCounts(obj) ########################
    df_EXP = np.log(obj.copy()+1)
    df_EXP['logCounts'] = np.log([smpl_counts[s.split('_')[0]] for s in df_EXP.index.values]) ########################
    regressors = {}
    for gene in tqdm(genes, total = len(genes)):
        regressor = SKLEARN_REGRESSOR_FACTORY[regressor_type](**regressor_kwargs)
        predictor_TF = getPredictorTFs(gene) #list(set(eRegulon_metadata.query("Gene == @gene")["TF"]))
        #remove gene itself as predictor
        if gene in predictor_TF:
            predictor_TF.remove(gene)
        if len(predictor_TF) == 0:
            continue
        predictor_TF_exp_v = df_EXP[predictor_TF+['logCounts']].to_numpy() ########################
        if len(predictor_TF_exp_v.shape) == 1:
            predictor_TF_exp_v = predictor_TF_exp_v.reshape(-1, 1)
        predictand_target_gene_exp_v = df_EXP[gene].to_numpy()
        regressor.fit(predictor_TF_exp_v, predictand_target_gene_exp_v) #ANCHOR - This is where i would add the log(nuclei_counts) as a predictor
        regressors[gene] = [*predictor_TF, regressor]
    #save regressors
    with open(save, "wb") as f:
        pkl.dump(regressors, f)
    return regressors

#map perturbed expression to original microglia PCA/UMAP space
# https://www.notion.so/loganlabbook/Sex-Differences-d379f0f2cd624cc0a505c12895158934?pvs=4#88f6d492b4f34af5a64a1fbf34b07317
# https://github.com/satijalab/seurat/blob/ca4c48b6300d7c10857e2d59a8ddee2858c7e2fc/R/dimensional_reduction.R#L282
# m = ProjectDim(m,reduction = 'umap')

print('loading data...')
obj = pd.read_csv(snakemake.input['cellStateMat'],index_col=0)
#filter down to genes in genesToTFs['mic0'].keys()
#obj = obj.loc[genesToTFs['mic0'].keys(),:]
obj = obj.T

print('data loaded')

reg_file = snakemake.output['regressors']
print('training regressors')
# if os.path.isfile(reg_file):
#     with open(reg_file, "rb") as f:
#         regressors = pkl.load(f)
# else:
regressors = train_gene_expression_models( obj = obj, save = reg_file)

# import pickle as pkl
# import numpy as np
# with open("/home/brasel/SingleCellProjects/MyProjects/pySCENIC_APOE/microglia/data/perturbation_mic01.pkl", "rb") as f:
#     perturbed_matrix = pkl.load(f)

#pca_loadings = pd.read_csv('microglia/data/PCA_feature_loadings.csv',index_col=0)
##len([g for g in pca_loadings.index.values if g in genes]) 1488 of the 3000 genes are in the data loaded into pySCENIC
#pca_loadings = pca_loadings.loc[[g for g in pca_loadings.index.values if g in obj.columns.values],:]

# #create a scatter plot of the umap object using columns "UMAP_1" and "UMAP_2" as the x and y coordinates
# import matplotlib.pyplot as plt
# plt.scatter(umap.iloc[:,0], umap.iloc[:,1], s=1)
# plt.savefig('microglia/plots/umap.pdf')
# plt.clf()
# plt.scatter(umap_pert.iloc[:,0], umap_pert.iloc[:,1], s=1)
# plt.savefig('microglia/plots/umap_pert.pdf')
# plt.clf()
# 
# 
# 
# 
# #Test how to replicated UMAP
# library(Seurat)
# library(Matrix)
# m = readRDS('/home/brasel/SingleCellProjects/dataObjects/microglia.rds')
# redeuc = m[['umap']]
# data.use = GetAssayData(m,slot='data')
# cell.embeddings <- Embeddings(object = redeuc)
# new.feature.loadings.full <- data.use %*% cell.embeddings
# 
# filt = rownames(GetAssayData(m,slot = 'scale.data'))
# test = data.frame(t(data.use[filt,]) %*% new.feature.loadings.full[filt,])
# test = data.frame(t(data.use) %*% new.feature.loadings.full)
# test$cluster = Idents(m)
# 
# library(ggplot2)
# p = ggplot(test,aes(x=UMAP_1,y=UMAP_2,color=cluster)) + geom_point(size=0.1)
# 
# 
# feature_loadings = m@reductions$pca@feature.loadings
# test = t(data.use[filt,]) %*% feature_loadings[filt,]
# m[['npca']] = m[['pca']]
# test = rbind(m[['npca']]@cell.embeddings,test)
# m[['npca']]@cell.embeddings = test[3001:6000,]
# 
# 
# obj = read.csv('/home/yanboyu/1PyScenicProject/discovery/DataObject/mic4.csv')
# filt = rownames(obj)
# test = data.frame(t(data.use[filt,]) %*% new.feature.loadings.full[filt,])
# #test = data.frame(t(data.use) %*% new.feature.loadings.full)
# test$cluster = Idents(m)
# 
# library(ggplot2)
# p = ggplot(test,aes(x=UMAP_1,y=UMAP_2,color=cluster)) + geom_point(size=0.1)
# 
# 
# ## plot umap coordinates from python
# umap = read.csv('microglia/data/umap.csv',row.names=1)
# umap_pert = read.csv('microglia/data/umap_pert.csv',row.names=1)
# 
# umap$cluster = Idents(m)[rownames(umap)]
# umap_pert$cluster = Idents(m)[rownames(umap_pert)]
# umap$pert = F
# umap_pert$pert = T
# 
# rownames(umap_pert) = paste0(rownames(umap_pert),'_pert')
# df = rbind(umap,umap_pert)
# 
# library(ggplot2)
# p = ggplot(df[df$cluster == 0,],aes(x=UMAP_1,y=UMAP_2,color=pert)) + geom_point(size=0.1)


#     #R
#     library(Seurat)
# library(Matrix)
# m = readRDS('/home/brasel/SingleCellProjects/dataObjects/microglia.rds')
# redeuc = m[['umap']]
# data.use = GetAssayData(m,slot='data')
# cell.embeddings <- Embeddings(object = redeuc)
# new.feature.loadings.full <- data.use %*% cell.embeddings
# filt = rownames(read.csv('/home/yanboyu/1PyScenicProject/discovery/DataObject/mic4.csv',row.names=1))
# emb = t(data.use[filt,]) %*% new.feature.loadings.full[filt,]
# write.csv(emb, file='/home/brasel/SingleCellProjects/MyProjects/pySCENIC_APOE/microglia/data/umap_embeddings_all_mic.csv')