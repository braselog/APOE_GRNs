import pandas as pd
import glob
import os
import re

def split_and_get_first_element(x):
    return x.split("(")[0]

def read_and_process_csv(file_path):
    df = pd.read_csv(file_path)
    return df['V1'].apply(split_and_get_first_element).tolist()

m0 = list(set(read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic0_newTF_TF_Gene_sum_above80.csv")).intersection(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic0_colonna_newTF_TF_Gene_sum_above80.csv")))

m1 = list(set(read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic1_newTF_TF_Gene_sum_above80.csv")).intersection(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic1_colonna_newTF_TF_Gene_sum_above80.csv")))

m3 = list(set(read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic3_newTF_TF_Gene_sum_above80.csv")).intersection(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic3_colonna_newTF_TF_Gene_sum_above80_20231031.csv")))

TFs = {"mic0": m0, "mic1": m1, "mic3": m3}

paths = ["/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic0_newTF_onehot", 
         "/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic1_newTF_onehot", 
         "/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic3_newTF_onehot"]
paths = dict(zip(["mic0", "mic1", "mic3"], paths))

#m = pd.read_csv("/home/brasel/SingleCellProjects/dataObjects/mic_01.csv",index_col=0)
m = pd.read_csv("/home/yanboyu/1PyScenicProject/discovery/DataObject/mic4.csv",index_col=0)
gene = m.index.values

byGene = {}
#for g in gene:
for clust in paths.keys():
    byGene[clust] = {}
    for g in gene:
        byGene[clust][g] = {}
    files = glob.glob(os.path.join(paths[clust], '*.csv'))
    tfs = sorted(TFs[clust])
    files = [f for f in files if any(tf in f for tf in tfs)]
    for f in files:
        df = pd.read_csv(f, index_col=0)
        for g in df.index:
            byGene[clust][g][re.split('_|\(',f)[3]] = df.loc[g, "sum"]

import pickle as pkl
with open(snakemake.output[[1]], "wb") as f:
    pkl.dump(byGene, f)

# #reload pickle data
# with open("/home/brasel/SingleCellProjects/MyProjects/pySCENIC_APOE/microglia/data/TF_hits_mic.pkl", "rb") as f:
#     test = pkl.load(f)