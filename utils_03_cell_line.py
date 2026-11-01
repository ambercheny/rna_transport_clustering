import pandas as pd
from pathlib import Path
import openpyxl
import re
import gzip
import re


def read_nucler_total_csv(final_gene_annotation, master_datasets, cell_line):

    annotation_lookup = final_gene_annotation[
        [
            "Original_Gene",
            "Current_Gene",
            "Ensembl_IDs",
            "Gene_type",
            "Main_analysis",
            "Mapping_method"
        ]
    ].copy()


    # ---------- M397 nuclear ----------
    m397_nuclear_annotated = (
        master_datasets[f"{cell_line}_nuclear"]
        .merge(
            annotation_lookup,
            left_on="Gene",
            right_on="Original_Gene",
            how="left"
        )
    )


    # ---------- M397 total ----------
    m397_total_annotated = (
        master_datasets[f"{cell_line}_total"]
        .merge(
            annotation_lookup,
            left_on="Gene",
            right_on="Original_Gene",
            how="left"
        )
    )


    print("M397 nuclear rows:", len(m397_nuclear_annotated))
    print("M397 total rows:  ", len(m397_total_annotated))

    print("\nNuclear annotation status:")
    print(f"m397_nuclear_annotated: shape {m397_nuclear_annotated.shape}, Gene.nunique() {m397_nuclear_annotated.Gene.nunique()}")
    print(m397_nuclear_annotated["Main_analysis"].value_counts(dropna=False))

    print("\nTotal RNA annotation status:")
    print(f"m397_total_annotated: shape {m397_total_annotated.shape}, Gene.nunique() {m397_total_annotated.Gene.nunique()}")
    print(m397_total_annotated["Main_analysis"].value_counts(dropna=False))

    print(f"Keep rows with Meta_analysis == 'KEEP'")

    m397_nuclear_pc = m397_nuclear_annotated[m397_nuclear_annotated["Main_analysis"] == "KEEP"].copy()
    m397_total_pc = m397_total_annotated[m397_total_annotated["Main_analysis"] == "KEEP"].copy()

    print("After only Meta_analysis = KEEP, M397 nuclear protein-coding:",
        f"{len(m397_nuclear_pc):,}")

    print("M397 total protein-coding:",
        f"{len(m397_total_pc):,}")

    print(f"m397_nuclear_pc: shape {m397_nuclear_pc.shape}, Current_Gene {m397_nuclear_pc.Current_Gene.nunique()}, Gene {m397_nuclear_pc.Gene.nunique()}")
    print(f"m397_total_pc: shape {m397_total_pc.shape}, Current_Gene {m397_total_pc.Current_Gene.nunique()}, Gene {m397_total_pc.Gene.nunique()}")
    return m397_nuclear_pc, m397_total_pc



def check_and_remove_duplicated_current_genes(m397_nuclear_pc, m397_total_pc, filename):
    # Identify collided current genes
    nuclear_collision_genes = set(
        m397_nuclear_pc.loc[
            m397_nuclear_pc["Current_Gene"].duplicated(keep=False),
            "Current_Gene"
        ]
    )

    total_collision_genes = set(
        m397_total_pc.loc[
            m397_total_pc["Current_Gene"].duplicated(keep=False),
            "Current_Gene"
        ]
    )

    print("Nuclear mapping-collision genes:",
        len(nuclear_collision_genes))

    print("Total mapping-collision genes:",
        len(total_collision_genes))


    # High-confidence one-to-one protein-coding tables
    m397_nuclear_pc_hc = m397_nuclear_pc[
        ~m397_nuclear_pc["Current_Gene"].isin(nuclear_collision_genes)
    ].copy().reset_index(drop=True)

    m397_total_pc_hc = m397_total_pc[
        ~m397_total_pc["Current_Gene"].isin(total_collision_genes)
    ].copy().reset_index(drop=True)


    print("\nHigh-confidence M397 nuclear PC genes:",
        f"{len(m397_nuclear_pc_hc):,}")

    print("High-confidence M397 total PC genes:",
        f"{len(m397_total_pc_hc):,}")

    print("\nDuplicate Current_Gene remaining:")
    print("Nuclear:",
        m397_nuclear_pc_hc["Current_Gene"].duplicated().sum())
    print("Total:",
        m397_total_pc_hc["Current_Gene"].duplicated().sum())


    # save collision genes as qc files
    collision_qc = pd.concat([
        m397_nuclear_pc[
            m397_nuclear_pc["Current_Gene"].isin(nuclear_collision_genes)
        ].assign(Dataset="M397_nuclear"),

        m397_total_pc[
            m397_total_pc["Current_Gene"].isin(total_collision_genes)
        ].assign(Dataset="M397_total")
    ])
    print(f"collision_qc: {collision_qc.shape}")
    
    # save collusion gene to qc file
    if not filename.exists():
        collision_qc.to_csv(filename,index=False)
    else:
        print(f"Use existing {filename}")
    return m397_nuclear_pc_hc, m397_total_pc_hc


def find_overlap(m397_nuclear_pc_hc, m397_total_pc_hc):
    nuclear_hc = set(m397_nuclear_pc_hc["Current_Gene"])
    total_hc = set(m397_total_pc_hc["Current_Gene"])

    shared_hc = nuclear_hc & total_hc
    nuclear_only_hc = nuclear_hc - total_hc
    total_only_hc = total_hc - nuclear_hc

    print("FINAL HIGH-CONFIDENCE OVERLAP")
    print("=" * 55)

    print("Nuclear PC genes:", f"{len(nuclear_hc):,}")
    print("Total PC genes:  ", f"{len(total_hc):,}")
    print("Shared genes:    ", f"{len(shared_hc):,}")
    print("Nuclear only:    ", f"{len(nuclear_only_hc):,}")
    print("Total only:      ", f"{len(total_only_hc):,}")

    print(
        "\n% nuclear genes represented in total:",
        f"{len(shared_hc) / len(nuclear_hc) * 100:.1f}%"
    )
