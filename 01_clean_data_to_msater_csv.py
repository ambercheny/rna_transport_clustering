from pathlib import Path

from utils_01_clean_data import inspect_excel, separate_excel_sheets, \
    check_gene_strings, save_duplicate_counts, standardize_gene_name, split_other_mln,\
    build_master_excel, check_master_integrity, save_master, build_project_gene_list

# QQQ we find those duplicated genes, but we still save the entire original data? yes

# Project paths
project_dir = Path.cwd().parent
print(f"project_dir: {project_dir}")
raw_dir = project_dir / "01_Raw_data"
clean_dir = project_dir / "02_Clean_data"

# Raw Excel files
nuclear_file = raw_dir / "Time series nuclear RNA expression.xlsx"
total_file = raw_dir / "Time series RNA-seq.xlsx"
print("Nuclear RNA file:", nuclear_file.exists())
print("Total RNA file:", total_file.exists())

# Read data
inspect_excel(nuclear_file, total_file)
datasets = separate_excel_sheets(nuclear_file, total_file)

# QC
duplicate_datasets = check_gene_strings(datasets)
save_duplicate_counts(duplicate_datasets, project_dir)

### prepare clean data
# standardize gene names
clean_datasets = standardize_gene_name(datasets)
# split MLN datasets
mln_split = split_other_mln(clean_datasets)
# create and integrity check master file
master_datasets = build_master_excel(clean_datasets, mln_split)
check_master_integrity(master_datasets)
# save clean data
save_master(clean_dir, master_datasets)

qc_dir = project_dir / "05_QC"
build_project_gene_list(master_datasets, qc_dir)