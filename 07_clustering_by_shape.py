import numpy as np
import pandas as pd
import pickle
from scipy.spatial.distance import pdist, squareform

from utils_clustering import get_trajectory_set_with_threshold, plot_centered_trajectory, \
    prepare_clustering_matrix, run_clustering
from unit_normalize_gene_vector import unit_normalize_gene_vector


# Unit-Normalize Each Gene

# ============================================================
# Analysis 2 — A2.5
# Unit-normalize centered log1p trajectories
# Focus clustering more strongly on trajectory shape
# ============================================================




def main():

    cell_line = "m397_nuclear_pc_hc"

    X = pd.read_csv(f"results/clustering_matrix_X_{cell_line}.csv", index_col=0)
    X_unit = unit_normalize_gene_vector(X)
    time_cols = ["0h", "12h", "36h", "72h"]
    k_evaluation, cluster_summaries = run_clustering(X_unit, time_cols, k_min = 2,k_max = 11,\
        figstr = f"_{cell_line}_unit")

if __name__ == "__main__":
    main()
    