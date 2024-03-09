import os
import glob
import pickle
import pandas as pd

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

'''
https://resources.aertslab.org/cistarget/
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
'''

def main(snakemake):
    DATABASES_GLOB = 'data/hg38_*genes_vs_motifs.rankings.feather'
    MOTIF_ANNOTATIONS_FNAME = 'data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
    MM_TFS_FNAME = 'data/allTFs_hg38.txt'
    MM_TFS_FNAME2 = 'data/allHumanTFs.txt'

    SC_EXP_FNAME = str(snakemake.input['expMat'])
    REGULONS_FNAME = snakemake.output['regulon_file']# "/home/yanboyu/1PyScenicProject/discovery/ResultedData/dis_astro0_oldTF/test_dis_astro0_regulonsL"+str(x)+".p"
    MOTIFS_FNAME = snakemake.output['motif_file']#"/home/yanboyu/1PyScenicProject/discovery/ResultedData/dis_astro0_oldTF/test_dis_astro0_motifsL"+str(x)+".csv"
    AUC_FNAME = snakemake.output['auc_file']

    if snakemake.params['run'] == 'first':
        # read in TFs
        tf_names = load_tf_names(MM_TFS_FNAME)
        tf_names2 = load_tf_names(MM_TFS_FNAME2)
        tmp = set(tf_names)
        tmp.update(tf_names2)
        tf_names = list(tmp)
        del tmp
        del tf_names2
    else:
        tf_names = load_tf_names(str(snakemake.input['stateTFs']))

    # read in single cell data
    ex_matrix = pd.read_csv(SC_EXP_FNAME, header=0, index_col=0).T
    ex_matrix.shape


    # load ranking databases
    db_fnames = glob.glob(DATABASES_GLOB)


    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]

    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
    dbs

    # infer co-expression modules
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True, seed=int(snakemake.params['seed']))

    # create modules from adjacencies
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)

    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)

    auc_mtx = aucell(ex_matrix, regulons, num_workers=4,seed=int(snakemake.params['seed']))
    auc_mtx.to_csv(AUC_FNAME)


if __name__ == '__main__':
    main(snakemake)

# import os
# import glob
# import pickle
# import pandas as pd
# import argparse
# 
# from dask.diagnostics import ProgressBar
# 
# from arboreto.utils import load_tf_names
# from arboreto.algo import grnboost2
# 
# from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
# from pyscenic.utils import modules_from_adjacencies, load_motifs
# from pyscenic.prune import prune2df, df2regulons
# from pyscenic.aucell import aucell
# 
# 
# def main(args):
#     DATABASES_GLOB = 'data/hg38_*genes_vs_motifs.rankings.feather'
#     MOTIF_ANNOTATIONS_FNAME = 'data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
#     MM_TFS_FNAME = 'data/allTFs_hg38.txt'
#     MM_TFS_FNAME2 = 'data/allHumanTFs.txt'
# 
#     SC_EXP_FNAME = args.expMat
#     REGULONS_FNAME = args.regulon_file
#     MOTIFS_FNAME = args.motif_file
#     AUC_FNAME = args.auc_file
# 
#     if args.run == 'first':
#         # read in TFs
#         tf_names = load_tf_names(MM_TFS_FNAME)
#         tf_names2 = load_tf_names(MM_TFS_FNAME2)
#         tmp = set(tf_names)
#         tmp.update(tf_names2)
#         tf_names = list(tmp)
#         del tmp
#         del tf_names2
#     else:
#         tf_names = load_tf_names(args.stateTFs)
# 
#     # read in single cell data
#     ex_matrix = pd.read_csv(SC_EXP_FNAME, header=0, index_col=0).T
#     ex_matrix.shape
# 
#     # load ranking databases
#     db_fnames = glob.glob(DATABASES_GLOB)
#     def name(fname):
#         return os.path.basename(fname).split(".")[0]
#     dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
# 
#     # co-expression modules
#     adjacencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names,seed=int(args.seed), verbose=True)
# 
#     # modules from adjacencies
#     modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
# 
#     # prune modules for targets
#     with ProgressBar():
#         df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
#     df.to_csv(MOTIFS_FNAME)
# 
#     # create regulons
#     regulons = df2regulons(df)
#     with open(REGULONS_FNAME, 'wb') as f:
#         pickle.dump(regulons, f)
# 
#     # calculate AUC matrix
#     auc_mtx = aucell(ex_matrix, regulons, num_workers=4,seed=int(args.seed))
#     auc_mtx.to_csv(AUC_FNAME)
# 
# if __name__ == "__main__":
#     # Argument parsing
#     parser = argparse.ArgumentParser(description='Run pySCENIC.')
#     parser.add_argument('--expMat', type=str, help='Expression Matrix')
#     parser.add_argument('--seed', type=str, help='Seed')
#     parser.add_argument('--run', type=str, help='Run')
#     parser.add_argument('--regulon_file', type=str, help='Regulon file')
#     parser.add_argument('--motif_file', type=str, help='Motif file')
#     parser.add_argument('--auc_file', type=str, help='AUC file')
#     parser.add_argument('--stateTFs', type=str, help='known cell state TFs')
# 
#     args = parser.parse_args()
#     main(args)