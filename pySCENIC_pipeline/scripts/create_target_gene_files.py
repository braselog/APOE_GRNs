# conda pySCENIC
import pickle
import pandas as pd

with open('snakemake_create_target_gene_files.pkl', 'wb') as f:
    pickle.dump(snakemake, f)

TFs_discovery = pd.read_csv(snakemake.input['TFs_discovery'], index_col=0)
TFs_replication = pd.read_csv(snakemake.input['TFs_replication'], index_col=0)
#TFs_discovery = pd.read_csv("mic/output/pySCENIC_discovery/second/mic0_discovery_TFs.csv", index_col=0)
#TFs_replication = pd.read_csv("mic/output/pySCENIC_discovery/second/mic0_discovery_TFs.csv", index_col=0)
TFs = set(TFs_discovery[TFs_discovery['sum'] >= 80].index.values).intersection(set(TFs_replication[TFs_replication['sum'] >= 50].index.values))

auc_files = snakemake.input['auc_mtx_files']
regulon_files = snakemake.input['regulon_files']
#file_names = "mic0_discovery_auc_mtx_0.csv mic0_discovery_auc_mtx_1.csv mic0_discovery_auc_mtx_2.csv mic0_discovery_auc_mtx_3.csv"
#auc_files = ["mic/output/pySCENIC_discovery/second/"+f for f in file_names.split()]
#file_names = "mic0_discovery_regulons_0.pkl mic0_discovery_regulons_1.pkl mic0_discovery_regulons_2.pkl mic0_discovery_regulons_3.pkl"
#regulon_files = ["mic/output/pySCENIC_discovery/second/"+f for f in file_names.split()]

#initialize a dictionary with the TFs as keys to store the dataframes
df_dict = {tf: pd.DataFrame() for tf in TFs}

for f,r in zip(auc_files,regulon_files):
    auc_mtx_tfs = [s[:-3] for s in pd.read_csv(f, index_col=0).columns.values]
    with open(r, "rb") as rfile:
        regulons = pickle.load(rfile)
    common_tfs = set(TFs).intersection(auc_mtx_tfs)

    for tf in common_tfs:
        df_dict[tf] = pd.merge(df_dict[tf],pd.DataFrame({f: [1 for x in range(len(regulons[auc_mtx_tfs.index(tf)].genes))]}, index = regulons[auc_mtx_tfs.index(tf)].genes), left_index=True, right_index=True,how='outer')

#replace all the NaN values with 0 and then add a sum column that sums the values of all the rows
for tf in df_dict:
    df_dict[tf] = df_dict[tf].fillna(0)
    df_dict[tf]['sum'] = df_dict[tf].sum(axis=1)/TFs_discovery.loc[tf,'sum']*100/df_dict[tf].shape[1]*100 #normalize by the number of times the TF showed up and then convert to a percentage
    df_dict[tf] = df_dict[tf].sort_values(by='sum', ascending=False)
    df_dict[tf].to_csv(str(snakemake.params['outFilePath']) + tf + "_target_genes.csv")

with open(snakemake.output[0], "w") as f:
    f.write("completed successfully!")
