import pandas as pd
import numpy as np
import pickle as pkl

with open('snakemake_filter_TFs.pkl','wb') as f:
    pkl.dump(snakemake,f)

csv_list = [[s[:-3] for s in pd.read_csv(f).columns.values] for f in snakemake.input['auc_mtx_files']]
#file_names = "mic0_replication_auc_mtx_0.csv mic0_replication_auc_mtx_1.csv mic0_replication_auc_mtx_2.csv mic0_replication_auc_mtx_3.csv"
#file_list = file_names.split()
#csv_list = [[s[:-3] for s in pd.read_csv('mic/output/pySCENIC_replication/second/'+f,index_col=0).columns.values] for f in file_list]

df=pd.DataFrame(data=csv_list)
df.index = ['list'+str(i) for i in range(1, df.shape[0]+1)]
df=df.T
df=df.drop(index=0)

# Get unique ids
id = df.stack().unique()

# Convert TF by onehot
y = pd.DataFrame({col: np.isin(id, df[col]) for col in df}, index=id)

# Sum each row
y['sum'] = y.sum(axis=1)

#Normalize each sum by the number of times pySCENIC was run. Only needed to make the test run work
y['sum'] = y['sum']/len(snakemake.input['auc_mtx_files'])*100

# Sort by 'sum' in descending order
y = y.sort_values(by='sum', ascending=False)

# Remove rows with empty index
y = y.loc[y.index != '']

# Write to csv
y.to_csv(snakemake.output['stateTFs_onehot'])

threshold = 80
if snakemake.params['run'] == 'second':
    if snakemake.wildcards['study'] == 'replication':
        threshold = 50
pd.Series([tf.split('(')[0] for tf in y.loc[y['sum'] >= threshold].index.values]).to_csv(snakemake.output['stateTFs'],index=False,header=False)