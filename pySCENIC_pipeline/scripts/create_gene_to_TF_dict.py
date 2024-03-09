import pandas as pd
import glob
import os
import re
import pickle as pkl

#with open('snakemake_create_gene_to_TF_dict.pkl', "wb") as f:
#    pkl.dump(snakemake, f)


#def split_and_get_first_element(x):
#    return x.split("(")[0]

def read_and_process_csv(file_path):
    df = pd.read_csv(file_path,header=None)
    return df.iloc[:,0].tolist()

disc_lists = [x for x in snakemake.input['TF_list'] if 'discovery' in x]
rep_lists = [x for x in snakemake.input['TF_list'] if 'replication' in x]
cellType = snakemake.wildcards['cellType']
states = [re.split(f'{cellType}|_d',x)[3] for x in disc_lists]#{'mic': ['0']}#,'1','3'], 'astro': ['0','1']}

stateRepTFs = [list(set(read_and_process_csv(d)).intersection(
    read_and_process_csv(r))) for d,r in zip(disc_lists,rep_lists)]

TFs = dict(zip([cellType+s for s in states], stateRepTFs))

#FIXME This path is not correct
paths = [f"../{cellType}/output/pySCENIC_discovery/targetGenes_{cellType}{state}" for state in states]
paths = dict(zip([cellType+s for s in states], paths))

gene = pd.read_csv(f"../{cellType}/data/pySCENIC/{cellType}0_discovery.csv",index_col=0).index.values


byGene = {}
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
            byGene[clust][g][re.split('/|_|\(',f)[7]] = df.loc[g, "sum"]

import pickle as pkl
with open(snakemake.output[0], "wb") as f:
    pkl.dump(byGene, f)