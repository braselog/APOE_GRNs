import pandas as pd
import glob
import os
import re
import upsetplot
import matplotlib.pyplot as plt

def split_and_get_first_element(x):
    return x.split("(")[0]

def read_and_process_csv(file_path,study):
    df = pd.read_csv(file_path)
    if study == 'discovery':
        threshold = 80
        df['TF'] = df['V1']
        df['sum'] = df['V2']
    else:
        threshold = 50
    df = df[df['sum'] >= threshold]
    return df['TF'].apply(split_and_get_first_element).tolist()

# this has all the colonna TFs, even those below 80: /home/yanboyu/1PyScenicProject/duplication/ResultedData/mic0_colonna_newTF_100run_onehot_20231030.csv
m0 = list(set(read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic0_newTF_TF_Gene_sum_above80.csv",'discovery')).intersection(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic0_colonna_newTF_100run_onehot_20231030C.csv",'replication')))

m1 = list(set(read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic1_newTF_TF_Gene_sum_above80.csv",'discovery')).intersection(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic1_colonna_newTF_100run_onehot_20231030C.csv",'replication')))

m3 = list(set(read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic3_newTF_TF_Gene_sum_above80.csv",'discovery')).intersection(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic3_colonna_newTF_100run_onehot_20231031C.csv",'replication')))

a0 = list(set(read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/astro0_newTF_TF_Gene_sum_above80.csv",'discovery')).intersection(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/astro0_colonna_newTF_100run_onehot_20231023C.csv",'replication')))

a1 = list(set(read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/astro1_newTF_TF_Gene_sum_above80.csv",'discovery')).intersection(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/astro1_colonna_newTF_100run_onehot_20231023C.csv",'replication')))

TFs = {"mic0": m0, "mic1": m1, "mic3": m3, "astro0": a0, "astro1": a1}

upsetplot.plot(upsetplot.from_contents(TFs), orientation = 'horizontal',sort_by='-degree', show_counts='%d', element_size=30)
plt.savefig('additionalAnalyses/plots/pySCENIC_TF_intersection_upsetplot.pdf')
plt.close()

#################################################################
paths = ["/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic0_newTF_onehot", 
         "/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic1_newTF_onehot", 
         "/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic3_newTF_onehot",
         "/home/yanboyu/1PyScenicProject/discovery/ResultedData/astro0_newTF_onehot", 
         "/home/yanboyu/1PyScenicProject/discovery/ResultedData/astro1_newTF_onehot"]
#paths = ["../mic/data/pySCENIC/second/targetGenes_mic0","../mic/data/pySCENIC/second/targetGenes_mic1","../mic/data/pySCENIC/second/targetGenes_mic3"]

paths = dict(zip(["mic0", "mic1", "mic3","astro0","astro1"], paths))

m = pd.read_csv("/home/yanboyu/1PyScenicProject/discovery/DataObject/mic4.csv",index_col=0)
#m = pd.read_csv("../mic/data/mic1_discovery.csv",index_col=0)
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
with open(snakemake.output[0], "wb") as f:
    pkl.dump(byGene, f)

# this has all the colonna TFs, even those below 80: /home/yanboyu/1PyScenicProject/duplication/ResultedData/mic0_colonna_newTF_100run_onehot_20231030.csv

# m0 = list(set(read_and_process_csv("../astro/data/pySCENIC/second/astro0_TF_list.csv")).intersection(
#     read_and_process_csv("../astro/data/pySCENIC_colonna/second/astro0_colonna_TF_list.csv")))
# 
# m1 = list(set(read_and_process_csv("../astro/data/pySCENIC/second/astro1_TF_list.csv")).intersection(
#     read_and_process_csv("../astro/data/pySCENIC_colonna/second/astro1_colonna_TF_list.csv")))
TFs = {"astro0": m0, "astro1": m1}

#paths = ["../astro/data/pySCENIC/second/targetGenes_astro0","../astro/data/pySCENIC/second/targetGenes_astro1"]
paths = dict(zip(["astro0", "astro1"], paths))

m = pd.read_csv("/home/yanboyu/1PyScenicProject/discovery/DataObject/astro4.csv",index_col=0)