import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


from load_master_datasets import load_master_datasets


def run_descriptive(df, time_cols):


    print("=" * 65)
    print("M397 NUCLEAR RNA — BASIC QC")
    print("=" * 65)

    print("\nShape:")
    print(df.shape)

    print("\nData types:")
    print(df[time_cols].dtypes)


    # ------------------------------------------------------------
    # 2. Missing / zero / negative values
    # ------------------------------------------------------------

    print("\n" + "=" * 65)
    print("MISSING / ZERO / NEGATIVE VALUES")
    print("=" * 65)

    qc = pd.DataFrame(index=time_cols)

    qc["Missing"] = df[time_cols].isna().sum()
    qc["Zero"] = (df[time_cols] == 0).sum()
    qc["Zero_%"] = qc["Zero"] / len(df) * 100
    qc["Negative"] = (df[time_cols] < 0).sum()

    print(qc.round(2))


    # ------------------------------------------------------------
    # 3. Expression summary
    # ------------------------------------------------------------

    print("\n" + "=" * 65)
    print("EXPRESSION SUMMARY")
    print("=" * 65)

    summary = df[time_cols].describe(
        percentiles=[0.01, 0.05, 0.10, 0.25, 0.50,
                    0.75, 0.90, 0.95, 0.99]
    )

    print(summary.round(4))


    # ------------------------------------------------------------
    # 4. Maximum values
    # ------------------------------------------------------------

    print("\n" + "=" * 65)
    print("MAXIMUM VALUES")
    print("=" * 65)

    for col in time_cols:
        print(f"{col:>4}: {df[col].max():.4f}")

    return time_cols


def all_gene_correlation_between_timepoints(df, time_cols):
    ### across all genes
    print("\n" + "=" * 65)
    print("PEARSON CORRELATION")
    print("=" * 65)

    print(
        df[time_cols]
        .corr(method="pearson")
        .round(3)
    )

    print("\nSPEARMAN CORRELATION")
    print("=" * 65)

    print(
        df[time_cols]
        .corr(method="spearman")
        .round(3)
    )


def plot_raw_distribution(df, time_cols):
# ------------------------------------------------------------
# 6. Distribution plots
# IMPORTANT: plot original values first — no additional log
# ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(4, 3))

    for col in time_cols:
        values = df.loc[df[col].notna(), col]
        ax.hist(
            values,
            bins=100,
            histtype="step",
            linewidth=1.5,
            label=col
        )

    ax.set_xlabel("Expression value (as provided)")
    ax.set_ylabel("Number of genes (log10)")
    ax.set_yscale("log")
    ax.set_title("M397 Nuclear RNA — Expression Distributions")
    ax.legend()

    plt.tight_layout()
    plt.savefig("raw_distribution.png")


# ------------------------------------------------------------
# 7. Zoom into central 99% of distribution
# ------------------------------------------------------------

def plot_raw_dis_quantile(df, time_cols, quantile=0.99):
    upper = df[time_cols].stack().quantile(quantile)

    fig, ax = plt.subplots(figsize=(4, 3))

    for col in time_cols:
        values = df.loc[
            df[col].notna() & (df[col] <= upper),
            col
        ]

        ax.hist(
            values,
            bins=100,
            histtype="step",
            linewidth=1.5,
            label=col
        )

    ax.set_xlim(0, upper)
    ax.set_xlabel("Expression value (as provided)")
    ax.set_ylabel("Number of genes (log10)")
    ax.set_title("M397 Nuclear RNA — Distribution (≤99th percentile)")
    ax.set_yscale("log")
    ax.legend()

    plt.tight_layout()
    plt.savefig(f"raw_distribution_q{str(quantile)}.png")


def check_anchor_raw_values(df, time_cols,
                            anchors = ["FTH1","FTL","HSPB1","GPX1","GPX4","SREBF2"]):
# ------------------------------------------------------------
# Expression of biological anchors
# ------------------------------------------------------------
    anchor_table = (
        df[df["Gene"].isin(anchors)]
        [["Gene"] + time_cols]
        .set_index("Gene")
        .reindex(anchors)
    )

    print("\n" + "=" * 65)
    print("BIOLOGICAL ANCHORS — RAW VALUES")
    print("=" * 65)

    print(anchor_table)


def plot_sample_trajectory(master_datasets, cell_line=["M397_nuclear", "M397_total", "M229_total"], genes = ["FTH1", "FTL"]):

    for index in cell_line:
    # ---------- M397 nuclear ----------
        m397_nuc_anchor = (
            master_datasets[index]
            .loc[lambda x: x["Gene"].isin(genes),
                ]
            .set_index("Gene")
        )
        print(f"cell lin : {index}")
        print(m397_nuc_anchor)

        if "M397_total" in index:
            m397_nuc_anchor = m397_nuc_anchor[[
                        "Control", "D3", "D8", "D13", "D21",
                        "D29", "D33", "D38", "D59",
                        "DR4", "DR10", "DR15", "DR17", "DR30", "DR35"
                    ]]
            m397_nuc_anchor.T.plot(marker="o", figsize=(6, 3))
            plt.title(f"{index.replace('_',' ')} total RNA — FTH1 / FTL long-term + withdrawal")
            plt.xlabel("Time")
            plt.ylabel("Original expression value")
            plt.xticks(rotation=45)
            plt.tight_layout()
            # plt.show()
        else:
            m397_nuc_anchor.T.plot(marker="o", figsize=(4, 3))
            plt.title(f"{index.replace('_',' ')} — FTH1 / FTL early trajectory")
            plt.xlabel("Time")
            plt.ylabel("Original expression value")
            plt.tight_layout()
        plt.savefig(f"{index}_trajectory.png")



def main():
    project_dir = Path.cwd().parent
    print(f"project_dir: {project_dir}")
    clean_dir = project_dir / "02_Clean_data"
    master_datasets = load_master_datasets(clean_dir)

    gene_name = "M397_nuclear"
    time_cols = ["0h", "12h", "36h", "72h"]

    nuc = master_datasets[gene_name].copy()
    time_cols = run_descriptive(nuc, time_cols)

    all_gene_correlation_between_timepoints(nuc, time_cols)

    plot_raw_distribution(nuc, time_cols)
    plot_raw_dis_quantile(nuc, time_cols, quantile=0.99)
    check_anchor_raw_values(nuc, time_cols)
    plot_sample_trajectory(master_datasets, cell_line=["M397_nuclear", "M397_total", "M229_total"], genes = ["FTH1", "FTL"]) #"M397_nuclear", "M397_total", "M229_total"


if __name__ == "__main__":
    main()





