"""
use tslearn package 
"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from tslearn.clustering import TimeSeriesKMeans
from utils_08_ts_cluster import plot_clusters
from unit_normalize_gene_vector import unit_normalize_gene_vector


cluster_shared = True
if cluster_shared:
    cell_line = "M397_shared"
else:
    cell_line = "M397_nuclear"


# cell_line = "m397_nuclear_pc_hc"
if cluster_shared:
    time_cols = ["Control", "D3", "D8", "D13", "D21", "D29", "D33", "D38"]
else:
    time_cols = ["0h", "12h", "36h", "72h"]




# data = pd.read_csv("data/M397_total_master.csv")
data = pd.read_csv(f"results/clustering_matrix_Xm397_total_pc_hc.csv", index_col=0)
X_unit = unit_normalize_gene_vector(data)
print(f"X_unit: {X_unit}")
if cluster_shared:
    data_total = pd.read_csv(f"results/clustering_matrix_Xm397_nuclear_pc_hc.csv", index_col=0)
    # data_total["Current_Gene"] = data_total.index
    shared_genes = X_unit.index.intersection(data_total.index)
    X_unit = X_unit.loc[shared_genes].copy()
    print(f"X_unit with shared genes: {X_unit}")

def cluster(model, X_ts, X_scaled, 
            save_filename="tslearn_gene_clusters"):

    labels = model.fit_predict(X_ts)
    # Save results


    clustered = X_scaled.copy()
    clustered["cluster"] = labels
    # clustered["variance"] = gene_variance.loc[clustered.index]

    savepath = f"{save_filename}.csv"
    clustered.to_csv(savepath, index=True, index_label="Gene")
    print(savepath)
    return clustered


n_clusters = [3] #[2,3,4,5,6,7,8,9,10]
for n in n_clusters:
        
    model = TimeSeriesKMeans(
        n_clusters=n,
        metric="euclidean",      # use "euclidean" if timing is strictly aligned
        max_iter=10,
        n_init = 50,
        random_state=42
    )
    ### time series clustering
    clustered = cluster(model, X_unit, X_unit,
                save_filename=f"results/tslearn_gene_clusters_{cell_line}_k{n}")

    plot_clusters(time_cols, clustered, X_unit,
                    fig_filename=f"results/cluster_results_{cell_line}_k{n}")