import pandas as pd
import pickle

from utils_clustering import get_trajectory_set_with_threshold, plot_centered_trajectory, \
    prepare_clustering_matrix, run_clustering

### after choosing a threshold from 05.py

# ============================================================
# A5.6 — Compare trajectory representations
# Primary working gene set: max expression > 0.01
# ============================================================


def main():

    threshold = 0.01
    cell_line = "m397_nuclear_pc_hc"
    m397_nuclear_pc_hc_expression = pd.read_csv(f"results/{cell_line}_expression.csv")
    time_cols = ["0h", "12h", "36h", "72h"]
    # time_cols = ["Control", "D3", "D8", "D13", "D21", "D29", "D33", "D38"]

    ### use threshold to filter genes and centering
    traj_raw, traj_centered = get_trajectory_set_with_threshold(m397_nuclear_pc_hc_expression, time_cols, threshold=threshold)
    plot_centered_trajectory(traj_centered, time_cols)

    ### clustering
    X = prepare_clustering_matrix(traj_centered, time_cols, X_str=cell_line)
    k_evaluation, cluster_summaries = run_clustering(X, time_cols, k_min = 2,k_max = 11, figstr=cell_line)

    # save results
    k_evaluation.to_csv(f"results/k_evaluation_{cell_line}.csv", index=False)
    with open(f"results/cluster_summaries_{cell_line}.pkl", "wb") as f:
        pickle.dump(cluster_summaries, f)


if __name__ == "__main__":
    main()
