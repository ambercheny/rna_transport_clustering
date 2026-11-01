import matplotlib.pyplot as plt
import math


def plot_clusters(time_cols, clustered, X_scaled,
                  fig_filename="cluster_results"):
    time_labels = time_cols
    clusters = sorted(clustered["cluster"].unique())
    n_clusters = len(clusters)
    n_cols = 5
    n_rows = math.ceil(n_clusters / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * 3, n_rows * 2),
        squeeze=False
    )
    axes = axes.ravel()

    for ax, c in zip(axes, clusters):
        genes = clustered.index[clustered["cluster"] == c]
        values = X_scaled.loc[genes, time_cols]

        # plot individual genes lightly
        for _, row in values.iterrows():
            ax.plot(time_labels, row, color="gray", alpha=0.08)

        # plot cluster mean
        ax.plot(
            time_labels,
            values.mean(axis=0),
            color="red",
            linewidth=2,
            label=f"Cluster {c}"
        )

        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_title(f"Cluster {c} n={len(genes)}")
        ax.tick_params(axis="x", rotation=45)
        ax.set_ylabel("Z-scored expression")

    for ax in axes[n_clusters:]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(f"{fig_filename}.pdf")
    plt.close(fig)


