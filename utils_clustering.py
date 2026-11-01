import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score
)
import pickle
from sklearn.metrics import adjusted_rand_score




def get_trajectory_set_with_threshold(df, time_cols, threshold=0.01):
    primary_mask = df["Max_expression"] > threshold

    traj_raw = (
        df.loc[
            primary_mask,
            ["Gene", "Current_Gene"] +time_cols
        ]
        .copy()
    )

    traj_log1p = np.log1p(traj_raw[time_cols])

    traj_z = traj_log1p.sub(
        traj_log1p.mean(axis=1),
        axis=0
    ).div(
        traj_log1p.std(axis=1),
        axis=0
    )

    print("Primary working gene set:", len(traj_raw))

    print("\nRaw:")
    print(traj_raw.describe())

    print("\nlog1p:")
    print(traj_log1p.describe())

    print("\nGene-wise z-score of log1p:")
    print(traj_z.describe())



    ### sanity check log center
    # gene-wise centering
    traj_centered = traj_log1p.sub(
        traj_log1p.mean(axis=1),
        axis=0
    )
    # Put gene names back
    traj_centered.insert(0, "Current_Gene", traj_raw["Current_Gene"].values)
    traj_centered.insert(0, "Gene", traj_raw["Gene"].values)
    print("Primary genes:", len(traj_centered))
    # Check that each gene has mean ~ 0 after centering
    row_means = traj_centered[time_cols].mean(axis=1)
    print("Largest absolute row mean after centering:",row_means.abs().max()) 
        # If centering worked correctly, every gene's mean must be identically/near zero.


    return traj_raw, traj_centered

def plot_centered_trajectory(
    traj_centered, time_cols,
    check_genes = ["FTH1","FTL","SPRY4","ETV4",],
    ):

    check_centered = (
        traj_centered[traj_centered["Current_Gene"].isin(check_genes)]
        .set_index("Current_Gene")[time_cols]
    )

    print("Centered log1p trajectories:")
    print(check_centered.round(4))


    check_centered.T.plot(
        marker="o",
        figsize=(8, 5)
    )

    plt.axhline(0, linewidth=1)
    plt.title("Centered log1p trajectories — sanity check")
    plt.xlabel("Time")
    plt.ylabel("Centered log1p expression")
    plt.tight_layout()
    plt.show()



def prepare_clustering_matrix(traj_centered, time_cols, X_str=""):
    # ============================================================
    # Analysis 2 — A2.1
    # Prepare centered log1p trajectory matrix
    # ============================================================

    X = (
        traj_centered
        .set_index("Current_Gene")[time_cols]
        .copy()
    )

    print("Clustering matrix shape:", X.shape)

    print("\nMissing values:", X.isna().sum().sum())
    print("Duplicated genes:", X.index.duplicated().sum())

    print(X.head())
    print(f"X type: {type(X)}, shape: {X.shape}")
    X.to_csv(f"results/clustering_matrix_X{X_str}.csv", index=True)
    return X


def plot_time_series_centroids(centroids, k, figstr=""):
    # Plot centroids
    centroids.T.plot(
        marker="o",
        figsize=(4, 3)
    )

    plt.axhline(0, linewidth=1)
    plt.title(f"K-means centroids — k={k}")
    plt.xlabel("Time")
    plt.ylabel("Centered log1p expression")
    plt.tight_layout()
    plt.savefig(f"results/cluster{k}{figstr}.png")


def run_clustering(X, time_cols, k_min = 2,k_max = 11, figstr=""):
    # ============================================================
    # Analysis 2 — A2.2
    # Evaluate candidate numbers of K-means clusters
    # ============================================================
    k_results = []
    cluster_summaries = {}

    for k in range(k_min, k_max):

        model = KMeans(n_clusters=k,random_state=42,n_init=50)

        labels = model.fit_predict(X)

        sizes = pd.Series(labels).value_counts().sort_index()
        centroids = pd.DataFrame(
                model.cluster_centers_,
                columns=time_cols
            )

        k_results.append({
            "k": k,
            "Silhouette": silhouette_score(X, labels),
            "Calinski_Harabasz": calinski_harabasz_score(X, labels),
            "Davies_Bouldin": davies_bouldin_score(X, labels),
            "Inertia": model.inertia_,
            "sizes":sizes,
        })
        cluster_summaries[k] = {
        "labels": labels,
        "sizes": sizes,
        "centroids": centroids
        }
        plot_time_series_centroids(centroids, k, figstr)

    k_evaluation = pd.DataFrame(k_results)

    print(k_evaluation.round(4))
    return k_evaluation, cluster_summaries


def run_multiple_clustering(X_unit):
    seeds = [1, 2, 3, 4, 5, 10, 20, 42, 100, 200]

    reference_model = KMeans(
        n_clusters=3,
        random_state=42,
        n_init=50
    )

    reference_labels = reference_model.fit_predict(X_unit)

    stability_results = []

    for seed in seeds:
        model = KMeans(
            n_clusters=3,
            random_state=seed,
            n_init=50
        )

        labels = model.fit_predict(X_unit)

        ari = adjusted_rand_score(
            reference_labels,
            labels
        )

        stability_results.append({
            "seed": seed,
            "Adjusted_Rand_Index": ari
        })

    stability_k = pd.DataFrame(stability_results)

    print(stability_k.round(4))

    return stability_k
