import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import numpy as np



def build_threshold_table(raw_pc):

    thresholds = [0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

    results_pc = []

    for threshold in thresholds:
        above = raw_pc > threshold

        results_pc.append({
            "threshold": threshold,
            ">=1_timepoint": (above.sum(axis=1) >= 1).sum(),
            ">=2_timepoints": (above.sum(axis=1) >= 2).sum(),
            ">=3_timepoints": (above.sum(axis=1) >= 3).sum(),
            "all_4_timepoints": (above.sum(axis=1) == 4).sum()
        })

    filter_sensitivity_pc = pd.DataFrame(results_pc)
    print(filter_sensitivity_pc)
    return filter_sensitivity_pc


def identify_peak_at_threshold(raw_pc, m397_nuclear_pc_hc, time_cols, threshold = 0.01):
    

    above = raw_pc > threshold
    n_above = above.sum(axis=1)

    single_tp = m397_nuclear_pc_hc.loc[
        n_above == 1,
        ["Gene", "Current_Gene"] + time_cols
    ].copy()

    # For each gene, identify peak time point (the time expression at peak) and peak value
    single_tp["Peak_time"] = single_tp[time_cols].idxmax(axis=1)
    single_tp["Peak_value"] = single_tp[time_cols].max(axis=1)

    # Second-highest value
    sorted_vals = np.sort(single_tp[time_cols].values, axis=1)
    single_tp["Second_highest"] = sorted_vals[:, -2]

    # Peak / second-highest ratio
    # this represents how high the peak compared to the second highest
    single_tp["Peak_to_second_ratio"] = (
        single_tp["Peak_value"] /
        single_tp["Second_highest"].replace(0, np.nan)
    )

    print("Single-timepoint-only genes:", len(single_tp))

    print("\nPeak time distribution:")
    print(single_tp["Peak_time"].value_counts())

    print("\nPeak-value summary:")
    print(single_tp["Peak_value"].describe())

    print("\nPeak-to-second-highest ratio summary:")
    print(single_tp["Peak_to_second_ratio"].describe())

    print("\nTop 20 strongest single-timepoint peaks:")
    print(
        single_tp
        .sort_values("Peak_to_second_ratio", ascending=False)
        .head(20)
    )


def max_expression(df, time_cols):
    # ============================================================
    # A5.4 — Max-expression and temporal dynamic range
    # High-confidence protein-coding genes
    # ============================================================
    # df = m397_nuclear_pc_hc
    df_expression = df[["Gene", "Current_Gene"]+time_cols].copy()

    # Maximum expression across the four time points
    df_expression["Max_expression"] = df_expression[time_cols].max(axis=1)
    # Minimum expression across the four time points
    df_expression["Min_expression"] = df_expression[time_cols].min(axis=1)
    print(f"df_expression")
    print(df_expression.head())

    # Absolute temporal dynamic range
    df_expression["Dynamic_range_raw"] = (df_expression["Max_expression"] - df_expression["Min_expression"])

    # Same idea after log1p transformation
    log_values = np.log1p(df_expression[time_cols])
    print("log_values")
    print(log_values.head())
    df_expression["Dynamic_range_log1p"] = (log_values.max(axis=1) - log_values.min(axis=1))

    print("Max-expression summary:")
    print(df_expression["Max_expression"].describe())

    print("\nRaw temporal dynamic-range summary:")
    print(df_expression["Dynamic_range_raw"].describe())

    print("\nLog1p temporal dynamic-range summary:")
    print(df_expression["Dynamic_range_log1p"].describe())

    return df_expression

def select_cutoff(df,
                  candidate_cutoffs = [0.005, 0.01, 0.02, 0.05]):
    # ============================================================
    # A5.5 — Evaluate candidate max-expression cutoffs
    # The purpose is to denoise, reducing the number of genes to analyze
    # ============================================================

    

    results = []

    for cutoff in candidate_cutoffs:

        keep = df["Max_expression"] > cutoff
        excluded = df.loc[~keep]

        results.append({
            "Max_cutoff": cutoff,
            "Genes_kept": keep.sum(),
            "Percent_kept": keep.mean() * 100,
            "Genes_excluded": (~keep).sum(),

            # check these dynamic range to answer, "Among the genes we are throwing away with Max_expression <= 0.01, did any of them have a meaningful fold change?"
            "Excluded_DR_median":
                excluded["Dynamic_range_log1p"].median(),

            "Excluded_DR_95th":
                excluded["Dynamic_range_log1p"].quantile(0.95),

            "Excluded_DR_max":
                excluded["Dynamic_range_log1p"].max()
        })

    cutoff_evaluation = pd.DataFrame(results)

    print(cutoff_evaluation.round(4))


def main():

    # define the preprocessed cell line to use
    # cell_line = "m397_nuclear_pc_hc"
    cell_line = "m397_total_pc_hc"
    m397_nuclear_pc_hc = pd.read_csv(f"results/{cell_line}.csv")
    # m397_total_pc_hc  = pd.read_csv(f"results/m397_total_pc_hc.csv")
    # time_cols = ["0h", "12h", "36h", "72h"]
    time_cols = ["Control", "D3", "D8", "D13", "D21", "D29", "D33", "D38"]
    raw_pc = m397_nuclear_pc_hc[time_cols].copy()
    filter_sensitivity_pc = build_threshold_table(raw_pc)
    identify_peak_at_threshold(raw_pc, m397_nuclear_pc_hc, time_cols, threshold = 0.01)
    m397_nuclear_pc_hc_expression = max_expression(m397_nuclear_pc_hc, time_cols)
    m397_nuclear_pc_hc_expression.to_csv(f"results/{cell_line}_expression.csv", index=False)
    select_cutoff(m397_nuclear_pc_hc_expression)


if __name__ == "__main__":
    main()


