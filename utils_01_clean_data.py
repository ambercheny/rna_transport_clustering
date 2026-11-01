import pandas as pd
from pathlib import Path
import openpyxl
import re

raw_gene_columns = {
    "M397 nuclear": "gene",
    "M397 total": "Gene",
    "M229 total": "Gene",
    "M263 total": "Gene",
    "Other MLN": "Genes",
    "HCC827 total": "Gene",
    "HT29 total": "Gene",
    "SKBR3 total": "Gene"
}

def inspect_excel(nuclear_file, total_file):
    # Inspect Excel workbook structure
    nuclear_xls = pd.ExcelFile(nuclear_file)
    total_xls = pd.ExcelFile(total_file)
    print("=== Nuclear RNA sheets ===")
    print(nuclear_xls.sheet_names)
    print("\n=== Total RNA sheets ===")
    print(total_xls.sheet_names)

    # Inspect the structure of every sheet

    print("=== NUCLEAR RNA ===")

    for sheet in nuclear_xls.sheet_names:
        df = pd.read_excel(nuclear_file, sheet_name=sheet)

        print(f"\nSheet: {sheet}")
        print("Shape:", df.shape)
        print("Columns:")
        print(df.columns.tolist())

        print("\nFirst 3 rows:")
        print(df.head(3))


    print("\n\n=== TOTAL RNA ===")

    for sheet in total_xls.sheet_names: # total_xls.sheet_names:
        df = pd.read_excel(total_file, sheet_name=sheet)

        print(f"\n{'='*70}")
        print(f"Sheet: {sheet}")
        print("Shape:", df.shape)

        print("Columns:")
        print(df.columns.tolist())

        print("\nFirst 3 rows:")
        print(df.head(3))



def separate_excel_sheets(nuclear_file, total_file):
    # Read all datasets with the correct header rows

    m397_nuclear = pd.read_excel(
        nuclear_file,
        sheet_name="M397 NUMI"
    )

    m397_total = pd.read_excel(
        total_file,
        sheet_name="M397 FPKM",
        header=1
    )

    m229_total = pd.read_excel(
        total_file,
        sheet_name="M229 FPKM"
    )

    m263_total = pd.read_excel(
        total_file,
        sheet_name="M263 FPKM"
    )

    other_mln = pd.read_excel(
        total_file,
        sheet_name="Other MLN lines Log2(FKPM+0.1)"
    )

    hcc827_total = pd.read_excel(
        total_file,
        sheet_name="HCC827 FPKM",
        header=1
    )

    ht29_total = pd.read_excel(
        total_file,
        sheet_name="HT-29 FPKM",
        header=1
    )

    skbr3_total = pd.read_excel(
        total_file,
        sheet_name="SKBR3 FPKM",
        header=1
    )

    datasets = {
        "M397 nuclear": m397_nuclear,
        "M397 total": m397_total,
        "M229 total": m229_total,
        "M263 total": m263_total,
        "Other MLN": other_mln,
        "HCC827 total": hcc827_total,
        "HT29 total": ht29_total,
        "SKBR3 total": skbr3_total
    }

    for name, df in datasets.items():
        print(f"{name:15s} {df.shape}")
        print(df.columns.tolist())
        print()
    return datasets

def save_duplicate_counts(duplicate_datasets, project_dir):
    qc_dir = project_dir / "05_QC"
    qc_dir.mkdir(exist_ok=True)
    clean_names = [k.replace("total", "").strip() for k in duplicate_datasets.keys()]
    dataset_tag = "_".join(clean_names)
    output_file = qc_dir / f"Duplicated_gene_symbols_{dataset_tag}.csv"

    
    duplicate_tables = []

    for name, df in duplicate_datasets.items():
        
        # Count how many times each gene symbol appears
        counts = df["Gene"].value_counts()
        
        # Keep only duplicated genes
        dup_counts = counts[counts > 1]
        
        # Create table
        temp = dup_counts.reset_index()
        temp.columns = ["Gene", "Number_of_rows"]
        temp.insert(0, "Cell_line", name)
        
        duplicate_tables.append(temp)
    
    if not duplicate_tables: 
        print("No duplicates found")
        return None
    else:
        # Combine all three cell lines
        duplicate_gene_summary = pd.concat(
            duplicate_tables,
            ignore_index=True
        )

        # Display complete table
        print(duplicate_gene_summary)

        # Save CSV
        if not output_file.exists():
            duplicate_gene_summary.to_csv(output_file, index=False)
            print(f"Saved to:\n{output_file}")

def check_gene_strings(datasets, raw_gene_columns=raw_gene_columns):                      
    print("GENE IDENTIFIER QC")                                                           
    print("=" * 75)                                                                       
                                                                                            
    duplicate_check = {}                                                                  
    duplicated_gene_sheet = []                                                            
                                                                                            
    for name, df in datasets.items():                                                     
        gene_col = raw_gene_columns[name]                                                 
        genes = df[gene_col].astype("string").str.strip()                                 
                                                                                            
        missing = genes.isna().sum()                                                      
        duplicates = genes.dropna().duplicated().sum()                                    
                                                                                            
        # Dynamically record datasets that contain duplicates                             
        if duplicates > 0:                                                                
            duplicated_gene_sheet.append(name)                                            
            duplicate_check[name] = (df, gene_col)  # Stores (DataFrame, gene_column_name)
                                                                                            
        antisense = genes.str.contains(r"-AS\d+$", regex=True, na=False).sum()            
        intronic = genes.str.contains(r"-IT\d+$", regex=True, na=False).sum()             
                                                                                            
        print(f"\n{name}")                                                                
        print(f"  Total rows:          {len(df):,}")                                      
        print(f"  Missing gene names:  {missing:,}")                                      
        print(f"  Duplicate genes:     {duplicates:,}")                                   
        print(f"  -AS# genes:          {antisense:,}")                                    
        print(f"  -IT# genes:          {intronic:,}")                                     
                                                                                            
    print("\nDatasets with duplicate genes:")                                             
    print(duplicated_gene_sheet)                                                          

    # convert to dictionary
    duplicate_datasets = {
            name: df.rename(columns={gene_col: "Gene"})
            for name, (df, gene_col) in duplicate_check.items()
        }
                                                           
    return duplicate_datasets  



def standardize_gene_name(datasets):
    # ============================================================
    # STEP 1: Standardize gene column names and remove whitespace
    # ============================================================

    # Make copies so the originally loaded data remain unchanged
    clean_datasets = {}

    for name, df in datasets.items():
        
        temp = df.copy()
        
        # Find the gene column used in the original dataset
        original_gene_col = raw_gene_columns[name]
        
        # Standardize the column name to "Gene"
        temp = temp.rename(columns={original_gene_col: "Gene"})
        
        # Remove leading/trailing whitespace from gene symbols
        temp["Gene"] = temp["Gene"].astype("string").str.strip()
        
        clean_datasets[name] = temp


    # Quick check
    for name, df in clean_datasets.items():
        
        missing = df["Gene"].isna().sum()
        duplicates = df[df["Gene"].duplicated()].Gene.nunique()
        
        # Check remaining leading/trailing spaces
        whitespace = (
            df["Gene"] != df["Gene"].str.strip()
        ).sum()

        print(
            f"{name:15s} | "
            f"rows={len(df):,}  gene column={df.columns[0]} | "
            f"missing={missing:3d} | "
            f"duplicates={duplicates:3d} | "
            f"whitespace={whitespace:3d}")
    return clean_datasets

def split_other_mln(clean_datasets):
    other_mln_clean = clean_datasets["Other MLN"]
    mln_cell_lines = [
        "M381",
        "M233",
        "M262",
        "M249",
        "M255",
        "M409"
    ]
    print(f"mln_cell_lines: {mln_cell_lines}")

    mln_split = {}

    for cell_line in mln_cell_lines:
        
        # Columns belonging to this cell line
        selected_columns = [
            "Gene",
            f"{cell_line}.Control",
            f"{cell_line}.D3",
            f"{cell_line}.D21"
        ]
        
        temp = other_mln_clean[selected_columns].copy()
        
        # Standardize time-point names
        temp = temp.rename(columns={
            f"{cell_line}.Control": "Control",
            f"{cell_line}.D3": "D3",
            f"{cell_line}.D21": "D21"
        })
        
        mln_split[cell_line] = temp


    # Inspect results
    for cell_line, df in mln_split.items():
        print(
            f"{cell_line}: "
            f"shape={df.shape}, "
            f"columns={df.columns.tolist()}"
        )
    return mln_split

def build_master_excel(clean_datasets, mln_split):
    # ============================================================
    # STEP 3: Build standardized master datasets
    # ============================================================

    master_datasets = {}

    # ---------- M397 nuclear ----------
    master_datasets["M397_nuclear"] = clean_datasets["M397 nuclear"].copy()


    # ---------- M397 total ----------
    master_datasets["M397_total"] = clean_datasets["M397 total"].copy()


    # ---------- M229 ----------
    m229 = clean_datasets["M229 total"].copy()

    m229 = m229.rename(columns={
        "M229_DMSO": "Control",
        "M229_3d": "D3",
        "M229_21d": "D21",
        "M229_60d": "D60"
    })

    master_datasets["M229_total"] = m229


    # ---------- M263 ----------
    m263 = clean_datasets["M263 total"].copy()

    m263 = m263.rename(columns={
        "M263_DMSO": "Control",
        "M263_3d": "D3",
        "M263_21d": "D21",
        "M263_90d": "D90"
    })

    master_datasets["M263_total"] = m263


    # ---------- Six other melanoma lines ----------
    for cell_line, df in mln_split.items():
        master_datasets[f"{cell_line}_total"] = df.copy()


    # ---------- Other cancer models ----------
    master_datasets["HCC827_total"] = clean_datasets["HCC827 total"].copy()
    master_datasets["HT29_total"] = clean_datasets["HT29 total"].copy()
    master_datasets["SKBR3_total"] = clean_datasets["SKBR3 total"].copy()


    # Check all datasets
    print(f"Total datasets: {len(master_datasets)}\n")

    for name, df in master_datasets.items():
        print(f"{name:15s} {df.shape}")
        print(f"  {df.columns.tolist()}")
    return master_datasets

def check_master_integrity(master_datasets):
    # ============================================================
    # STEP 4: Integrity check before saving master datasets
    # ============================================================

    print("MASTER DATA INTEGRITY CHECK")
    print("=" * 70)

    all_pass = True

    for name, df in master_datasets.items():

        problems = []

        # 1. Gene column must exist
        if "Gene" not in df.columns:
            problems.append("Missing Gene column")

        # 2. Gene names should not be missing
        if "Gene" in df.columns:
            missing_genes = df["Gene"].isna().sum()
            if missing_genes > 0:
                problems.append(f"{missing_genes} missing gene names")

        # 3. Check expression columns for unexpected missing values
        numeric_cols = df.select_dtypes(include="number").columns
        missing_values = df[numeric_cols].isna().sum().sum()

        if missing_values > 0:
            problems.append(f"{missing_values} missing expression values")

        # Report
        if problems:
            all_pass = False
            print(f"FAIL  {name}")
            for p in problems:
                print(f"      - {p}")
        else:
            print(f"PASS  {name:15s}  rows={len(df):,}")


    print("\n" + "=" * 70)

    if all_pass:
        print("ALL DATASETS PASSED ✓")
    else:
        print("Some datasets need inspection before saving.")



def save_master(clean_dir, master_datasets):
    # ============================================================
    # STEP 5: Save standardized master datasets
    # ============================================================

    clean_dir.mkdir(exist_ok=True)

    for name, df in master_datasets.items():

        output_file = clean_dir / f"{name}_master.csv"
        if not output_file.exists():

            df.to_csv(
                output_file,
                index=False
            )

            print(f"Saved: {output_file.name}")
            print("\nDone!")
            print(f"Master datasets saved to:\n{clean_dir}")
        else:
            print(f"Not saved; use existing clean master")



def build_project_gene_list(master_datasets, qc_dir):
    # ============================================================
    # STEP 6: Build project-wide unique gene list
    # ============================================================

    all_genes = []

    for dataset_name, df in master_datasets.items():
        
        genes = df["Gene"].dropna().astype(str).str.strip()
        
        for gene in genes:
            all_genes.append({
                "Gene": gene,
                "Dataset": dataset_name
            })

    gene_presence = pd.DataFrame(all_genes)

    # Number of datasets in which each gene appears
    gene_summary = (
        gene_presence
        .groupby("Gene")
        .agg(
            Number_of_datasets=("Dataset", "nunique"),
            Datasets=("Dataset", lambda x: "; ".join(sorted(set(x))))
        )
        .reset_index()
        .sort_values(
            ["Number_of_datasets", "Gene"],
            ascending=[False, True]
        )
    )

    print("Total gene entries across datasets:", f"{len(gene_presence):,}")
    print("Unique gene symbols:", f"{len(gene_summary):,}")

    print(gene_summary.head(20))

    # save
    qc_dir.mkdir(exist_ok=True)
    
    gene_list_file = qc_dir / "Project_unique_gene_symbols.csv"
    if not gene_list_file.exists():
        gene_summary.to_csv(
            gene_list_file,
            index=False
        )
        print(f"Saved:\n{gene_list_file}")
    else:
        print(f"Not saved; use existing {gene_list_file}")


