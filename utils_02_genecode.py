import pandas as pd
from pathlib import Path
import openpyxl
import re
import gzip
import re
import pandas as pd
from collections import defaultdict


def read_gencode(gencode_file):

    gene_records = []

    with gzip.open(gencode_file, "rt") as f:
        for line in f:

            # Skip comments
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            # GTF should contain 9 fields
            if len(fields) != 9:
                continue

            chromosome = fields[0]
            feature = fields[2]
            attributes = fields[8]

            # We only need gene-level entries
            if feature != "gene":
                continue

            gene_id = re.search(r'gene_id "([^"]+)"', attributes)
            gene_name = re.search(r'gene_name "([^"]+)"', attributes)
            gene_type = re.search(r'gene_type "([^"]+)"', attributes)

            gene_records.append({
                "Gene": gene_name.group(1) if gene_name else None,
                "Ensembl_ID": gene_id.group(1) if gene_id else None,
                "Gene_type": gene_type.group(1) if gene_type else None,
                "Chromosome": chromosome
            })


    gene_annotation = pd.DataFrame(gene_records)

    print("\nGENCODE gene-level entries:", f"{len(gene_annotation):,}")
    print("Unique gene symbols:", f"{gene_annotation['Gene'].nunique():,}")

    print(gene_annotation.head(10))

    ### sanity check
    check_genes = [
        "FTH1",
        "FTL",
        "HSPB1",
        "A1BG",
        "A1BG-AS1",
        "SNORA11"
    ]
    print("Sanity check")
    print(  
        gene_annotation[
            gene_annotation["Gene"].isin(check_genes)
        ].sort_values("Gene")
    )

    ### check Genecode with multiple Gene
    # ============================================================
    # STEP 8: Map project genes to GENCODE v50
    # ============================================================

    # Keep the annotation fields we need
    annotation_map = gene_annotation[["Gene", "Ensembl_ID", "Gene_type", "Chromosome"]].copy()

    # Check whether one gene symbol maps to multiple GENCODE entries
    annotation_counts = (
        annotation_map
        .groupby("Gene")
        .size()
    )

    multi_mapping_genes = annotation_counts[annotation_counts > 1]

    # these are usually abundant non-coding RNA genes; for example, 5S_rRNA present in multiple chromosomes
    # In GENECODE, every physical location gets its own Ensembl ID
    print(
        "GENCODE gene symbols with multiple entries; one Gene name with multiple Ensemble ID",
        f"{len(multi_mapping_genes):,}"
    )


    return gene_annotation

def merge_project_with_genecode(gene_summary, annotation_map):

    # Merge project gene list with GENCODE
    project_annotation = gene_summary.merge(
        annotation_map,
        on="Gene",
        how="left"
    )

    print("\nProject unique genes:", f"{len(gene_summary):,}")

    matched = project_annotation["Gene_type"].notna().sum()
    unmatched = project_annotation["Gene_type"].isna().sum()
    print("Mapped rows:", f"{matched:,}")
    print("Unmapped rows:", f"{unmatched:,}")
    print(f"="*20)

    mapped_genes = project_annotation[project_annotation["Gene_type"].notna()]["Gene"].nunique()
    unmapped_genes = project_annotation[project_annotation["Gene_type"].isna()]["Gene"].nunique()
    print(f"Total unique genes:   {len(gene_summary):,}") 
    print(f"Mapped unique genes:   {mapped_genes:,}")                                                         
    print(f"Unmapped unique genes: {unmapped_genes:,}")

    return mapped_genes, unmapped_genes, project_annotation

def genecode_mapping_status(project_annotation):

    mapping_summary = (
        project_annotation
        .groupby("Gene", as_index=False)
        .agg(
            Number_of_datasets=("Number_of_datasets", "first"),
            Datasets=("Datasets", "first"),
            Number_of_GENCODE_matches=("Ensembl_ID", "count"),
            Gene_types=("Gene_type",
                        lambda x: "; ".join(sorted(set(x.dropna())))),
            Ensembl_IDs=("Ensembl_ID",
                        lambda x: "; ".join(sorted(set(x.dropna())))),
            Chromosomes=("Chromosome",
                        lambda x: "; ".join(sorted(set(x.dropna()))))
        )
    )

    mapping_summary["Mapping_status"] = "Matched"

    mapping_summary.loc[
        mapping_summary["Number_of_GENCODE_matches"] == 0,
        "Mapping_status"
    ] = "Unmatched"

    mapping_summary.loc[
        mapping_summary["Number_of_GENCODE_matches"] > 1,
        "Mapping_status"
    ] = "Multiple"


    print(mapping_summary.groupby(["Mapping_status"]).Gene.nunique())

    print(
        "\nTotal project genes:",
        f"{len(mapping_summary):,}"
    )
    return mapping_summary


def check_gencode_unmatched_genes(mapping_summary, master_datasets):
    # check unmatched genes
    unmatched_genes = mapping_summary[
        mapping_summary["Mapping_status"] == "Unmatched"
    ].copy()

    print("Total unmatched genes:", f"{len(unmatched_genes):,}")

    # How many unmatched genes appear in each dataset?
    print("\nUnmatched genes by dataset:")
    print("=" * 50)

    for dataset_name in sorted(master_datasets.keys()):
        
        genes_in_dataset = set(
            master_datasets[dataset_name]["Gene"]
            .dropna()
            .astype(str)
            .str.strip()
        )
        
        unmatched_in_dataset = set(unmatched_genes["Gene"]) & genes_in_dataset
        
        print(
            f"{dataset_name:15s}: "
            f"{len(unmatched_in_dataset):,} / {len(genes_in_dataset):,} "
            f"({len(unmatched_in_dataset)/len(genes_in_dataset)*100:.1f}%)"
        )
    return unmatched_genes
        
    
def read_hgnc(hgnc_file):
    if hgnc_file.exists():
        print("File size:", round(hgnc_file.stat().st_size / 1024 / 1024, 1), "MB")
    # ============================================================
    # STEP 12: Load and inspect HGNC complete gene set
    # ============================================================

    hgnc = pd.read_csv(hgnc_file,sep="\t", dtype=str)

    print("HGNC rows:", f"{len(hgnc):,}")
    print("HGNC columns:", len(hgnc.columns))

    print("\nColumn names:")
    print(hgnc.columns.tolist())

    # Show key HGNC annotation fields
    key_columns = [
        "hgnc_id",
        "symbol",
        "name",
        "locus_group",
        "locus_type",
        "prev_symbol",
        "alias_symbol",
        "ensembl_gene_id"
    ]

    existing_columns = [
        col for col in key_columns
        if col in hgnc.columns
    ]

    print(hgnc[existing_columns].head(10))

    return hgnc


def test_hgnc(hgnc, test_symbol = "AARS"):

    # Search approved/current symbol
    current_match = hgnc[
        hgnc["symbol"] == test_symbol
    ]

    # Search previous symbols
    previous_match = hgnc[
        hgnc["prev_symbol"]
        .fillna("")
        .str.split("|")
        .apply(lambda x: test_symbol in x)
    ]

    # Search aliases
    alias_match = hgnc[
        hgnc["alias_symbol"]
        .fillna("")
        .str.split("|")
        .apply(lambda x: test_symbol in x)
    ]

    print("Current symbol matches:", len(current_match))
    print("Previous symbol matches:", len(previous_match))
    print("Alias matches:", len(alias_match))

    print(
        previous_match[
            [
                "hgnc_id",
                "symbol",
                "name",
                "locus_group",
                "locus_type",
                "prev_symbol",
                "alias_symbol",
                "ensembl_gene_id"
            ]
        ]
    )


def build_hgnc_symbol_map(hgnc):

    previous_map = defaultdict(set)
    alias_map = defaultdict(set)

    for _, row in hgnc.iterrows():
        
        current_symbol = row["symbol"]
        
        # Previous symbols
        if pd.notna(row["prev_symbol"]):
            for old_symbol in row["prev_symbol"].split("|"):
                old_symbol = old_symbol.strip()
                if old_symbol:
                    previous_map[old_symbol].add(current_symbol)
        
        # Alias symbols
        if pd.notna(row["alias_symbol"]):
            for alias in row["alias_symbol"].split("|"):
                alias = alias.strip()
                if alias:
                    alias_map[alias].add(current_symbol)


    print("Unique previous symbols:", f"{len(previous_map):,}")
    print("Unique alias symbols:", f"{len(alias_map):,}")
    return previous_map, alias_map


def label_with_hgnc_rescue(hgnc, unmatched_genes,previous_map,alias_map):

    hgnc_current_symbols = set(hgnc["symbol"].dropna())

    rescue_records = []

    for gene in unmatched_genes["Gene"]:
        
        current_gene = None
        method = None
        
        # Current HGNC symbol
        if gene in hgnc_current_symbols:
            current_gene = gene
            method = "HGNC_current"
        
        # Previous HGNC symbol
        elif gene in previous_map:
            
            matches = previous_map[gene]
            
            if len(matches) == 1:
                current_gene = next(iter(matches))
                method = "HGNC_previous"
            else:
                method = "ambiguous_previous"
        
        # HGNC alias
        elif gene in alias_map:
            
            matches = alias_map[gene]
            
            if len(matches) == 1:
                current_gene = next(iter(matches))
                method = "HGNC_alias"
            else:
                method = "ambiguous_alias"
        
        else:
            method = "unmatched"
        
        
        rescue_records.append({
            "Original_Gene": gene,
            "Current_Gene": current_gene,
            "Rescue_method": method
        })


    rescue_df = pd.DataFrame(rescue_records)

    print("RESCUE RESULTS")
    print("=" * 50)
    print(rescue_df["Rescue_method"].value_counts())

    print("\nTotal:", f"{len(rescue_df):,}")

    return rescue_df

def validate_genecode_with_hgnc(rescue_df, gene_annotation):

    # Only genes with a unique rescued current symbol
    rescued = rescue_df[
        rescue_df["Current_Gene"].notna()
    ].copy()

    # Join rescued current symbols back to GENCODE
    rescued_validated = rescued.merge(
        gene_annotation,
        left_on="Current_Gene",
        right_on="Gene",
        how="left"
    )

    # Remove redundant GENCODE Gene column
    rescued_validated = rescued_validated.drop(columns=["Gene"])

    print("Unique rescued original genes:",
        f"{rescued['Original_Gene'].nunique():,}")

    print("Rescued genes with GENCODE annotation:",
        f"{rescued_validated.loc[rescued_validated['Gene_type'].notna(), 'Original_Gene'].nunique():,}")

    print("Rescued genes still without GENCODE annotation:",
        f"{rescued_validated.loc[rescued_validated['Gene_type'].isna(), 'Original_Gene'].nunique():,}")

    return rescued_validated

def get_rescued_biotypes_from_hgnc(rescued_validated):
    # Gene biotypes recovered through HGNC rescue

    rescued_biotypes = (
        rescued_validated[
            rescued_validated["Gene_type"].notna()
        ]
        .drop_duplicates(["Original_Gene", "Current_Gene", "Gene_type"])
        ["Gene_type"]
        .value_counts()
    )

    print("\nRecovered GENCODE gene types:")
    print(rescued_biotypes)
    return rescued_biotypes


def build_master_gene_annotation(gene_summary, gene_annotation, rescued_validated):

    project_gene_master = gene_summary[
        ["Gene", "Number_of_datasets", "Datasets"]
    ].copy()
    project_gene_master = project_gene_master.rename(
    columns={"Gene": "Original_Gene"}
)
### 1. direct GENECODE
    direct_annotation = gene_annotation.rename(columns={
        "Gene": "Original_Gene"
    })

    direct_annotation["Current_Gene"] = direct_annotation["Original_Gene"]
    direct_annotation["Mapping_method"] = "GENCODE_direct"

    direct_annotation = direct_annotation[
        [
            "Original_Gene",
            "Current_Gene",
            "Ensembl_ID",
            "Gene_type",
            "Chromosome",
            "Mapping_method"
        ]
    ]

    # Only unique direct mappings for now
    direct_counts = direct_annotation["Original_Gene"].value_counts()

    direct_unique = direct_annotation[
        direct_annotation["Original_Gene"].map(direct_counts) == 1
    ].copy()

### 2. HGNC rescue
    rescue_annotation = (
        rescued_validated[
            rescued_validated["Gene_type"].notna()
        ]
        [
            [
                "Original_Gene",
                "Current_Gene",
                "Ensembl_ID",
                "Gene_type",
                "Chromosome",
                "Rescue_method"
            ]
        ]
        .copy()
    )

    rescue_annotation = rescue_annotation.rename(
        columns={"Rescue_method": "Mapping_method"}
    )

    # Keep only one identical annotation record
    rescue_annotation = rescue_annotation.drop_duplicates()

### 3. concat both
    combined_annotation = pd.concat(
        [
            direct_unique,
            rescue_annotation
        ],
        ignore_index=True
    )



    project_gene_master = project_gene_master.merge(
        combined_annotation,
        on="Original_Gene",
        how="left"
    )

    

### 4. define Main_analysis

    project_gene_master["Main_analysis"] = "EXCLUDE"

    project_gene_master.loc[
        project_gene_master["Gene_type"] == "protein_coding",
        "Main_analysis"
    ] = "KEEP"

    # Anything without a resolved annotation remains unresolved,
    # NOT automatically classified as non-coding.
    project_gene_master.loc[
        project_gene_master["Gene_type"].isna(),
        "Main_analysis"
    ] = "UNRESOLVED"

    print("MAIN ANALYSIS STATUS")
    print("=" * 50)

    print(project_gene_master["Main_analysis"].value_counts())

    print("\nGENE TYPES")
    print("=" * 50)

    print(
        project_gene_master["Gene_type"]
        .value_counts(dropna=False)
        .head(20)
    )

    return project_gene_master


def qc_duplicated_gene_rows(project_gene_master, gene_summary):
    master_counts = (
        project_gene_master["Original_Gene"]
        .value_counts()
    )

    duplicated_master_genes = master_counts[
        master_counts > 1
    ]

    print("Original genes with >1 annotation row:",
        len(duplicated_master_genes))

    print("Extra rows:",
        len(project_gene_master) - len(gene_summary))

    print("\nNumber of annotation rows per affected gene:")
    print(
        duplicated_master_genes
        .value_counts()
        .sort_index()
    )

    ### print some genes
    # check the duplicated gene
    affected = project_gene_master[
        project_gene_master["Original_Gene"].isin(
            duplicated_master_genes.index
        )
    ].sort_values(
        ["Original_Gene", "Current_Gene", "Ensembl_ID"]
    )

    print(
        affected[
            [
                "Original_Gene",
                "Current_Gene",
                "Ensembl_ID",
                "Gene_type",
                "Chromosome",
                "Mapping_method",
                "Main_analysis"
            ]
        ].head(100)
    )

    print(f"Duplicated information per Original_Gene")
    affected = project_gene_master[
        project_gene_master["Original_Gene"].isin(
            duplicated_master_genes.index
        )
    ].copy()
    print(f"affected Original_Gene: {affected.Original_Gene.nunique()}, Current_Gene: {affected.Current_Gene.nunique()}")


def build_final_gene_annotation(project_gene_master, gene_summary):
    # ============================================================
    # Collapse annotation to exactly one row per original gene; 
    # where if there are >1 ensemble IDs, current genes, it gets stored in the corresponding columns; no metadata is excluded
    # ============================================================

    def combine_unique(series):
        values = sorted(set(series.dropna().astype(str)))
        return "; ".join(values) if values else pd.NA


    final_gene_annotation = (
        project_gene_master
        .groupby("Original_Gene", as_index=False)
        .agg(
            Number_of_datasets=("Number_of_datasets", "first"),
            Datasets=("Datasets", "first"),
            Current_Gene=("Current_Gene", combine_unique),
            Ensembl_IDs=("Ensembl_ID", combine_unique),
            Gene_type_details=("Gene_type", combine_unique),
            Chromosomes=("Chromosome", combine_unique),
            Mapping_method=("Mapping_method", combine_unique)
        )
    )

    # Number of distinct gene types
    final_gene_annotation["Number_of_gene_types"] = (
        final_gene_annotation["Gene_type_details"]
        .fillna("")
        .apply(
            lambda x: len(x.split("; ")) if x else 0
        )
    )

    # Final simplified gene type
    final_gene_annotation["Gene_type"] = final_gene_annotation["Gene_type_details"]

    final_gene_annotation.loc[
        final_gene_annotation["Number_of_gene_types"] > 1,
        "Gene_type"
    ] = "AMBIGUOUS"


    # Main-analysis status
    final_gene_annotation["Main_analysis"] = "EXCLUDE"

    final_gene_annotation.loc[
        final_gene_annotation["Gene_type"] == "protein_coding",
        "Main_analysis"
    ] = "KEEP"

    final_gene_annotation.loc[
        final_gene_annotation["Gene_type"].isna(),
        "Main_analysis"
    ] = "UNRESOLVED"

    final_gene_annotation.loc[
        final_gene_annotation["Gene_type"] == "AMBIGUOUS",
        "Main_analysis"
    ] = "AMBIGUOUS"

    print("FINAL ANNOTATION QC")
    print("=" * 50)

    print("Original project genes:",
        f"{len(gene_summary):,}")

    print("Final annotation rows:",
        f"{len(final_gene_annotation):,}")

    print("Duplicate Original_Gene:",
        final_gene_annotation["Original_Gene"].duplicated().sum())

    print("\nMain-analysis status:")
    print(final_gene_annotation["Main_analysis"].value_counts())

    print("\nAmbiguous genes:",
        (final_gene_annotation["Main_analysis"] == "AMBIGUOUS").sum())

    return final_gene_annotation



def save_final_gene_annotation(final_gene_annotation, output_file):
# ============================================================
# STEP 21: Save final project gene annotation master
# ============================================================

    if not output_file.exists():
        final_gene_annotation.to_csv(
            output_file,
            index=False
        )
        print("Saved:")
        print(output_file)
        print()
        print("Rows:", f"{len(final_gene_annotation):,}")
    else:
        print(f"No save; keep existing {output_file}")


def save_unresolved_genes(final_gene_annotation, qc_file):

    if not qc_file.exists():
        final_gene_annotation[
            final_gene_annotation["Main_analysis"].isin(
                ["UNRESOLVED", "AMBIGUOUS"]
            )
        ].to_csv(
            qc_file,
            index=False
        )
        print("QC file saved:")
        print(qc_file)
    else:
        print(f"No save; keep existing QC {qc_file}")
