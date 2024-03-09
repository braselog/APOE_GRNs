import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from scipy.io import mmread
import pickle as pkl
import matplotlib.pyplot as plt

# Load the GEM and UMAP data from file
geneFilt = pd.read_csv(snakemake.input['geneFilt'],header=None)
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
#gem.columns = gem.columns.astype(str)
train_gem, test_gem, train_umat, test_umat = train_test_split(gem, umat, test_size=0.2, stratify=cluster)
train_cluster = cluster.loc[train_gem.index]
train_gem, val_gem, train_umat, val_umat = train_test_split(train_gem, train_umat, test_size=0.1, stratify=train_cluster)

# Normalize the GEM data
scaler = StandardScaler()
train_gem = scaler.fit_transform(train_gem)
val_gem = scaler.transform(val_gem)
test_gem = scaler.transform(test_gem)

# Save scaler
with open(snakemake.output['scaler'], 'wb') as f:
    pkl.dump(scaler, f)

# Define the neural network architecture
model = tf.keras.models.Sequential([
    tf.keras.layers.Dense(128, activation='relu', input_dim=train_gem.shape[1]),
    tf.keras.layers.Dense(64, activation='relu'),
    tf.keras.layers.Dense(train_umat.shape[1])
])

# Compile the model
model.compile(optimizer='adam', loss='mse')

# Define early stopping
early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)

# Train the model
model.fit(train_gem, train_umat, validation_data=(val_gem, val_umat), epochs=1000, callbacks=[early_stop])

# Evaluate the model
test_loss = model.evaluate(test_gem, test_umat)
print('Test loss:', test_loss)

# Compare the predicted UMAT values to the actual UMAT values by plotting them
new_umat = model.predict(scaler.transform(gem))

# Save new_umat to csv file
new_umat = pd.DataFrame(new_umat, index=gem.index)
new_umat.to_csv(snakemake.output['UMAP_cell_embeddings'])

# Save model
model.save(snakemake.output['model'])