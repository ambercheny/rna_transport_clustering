"""Per cell line, do the preprocessing. e.g., remove Current_gene that has more than one legacy entity"""
"""Use M397 as example here"""
import pandas as pd
from pathlib import Path
import openpyxl
import re
import gzip
import re
project_dir = Path.cwd().parent
print(f"project_dir: {project_dir}")

from load_master_datasets import load_master_datasets
from utils_03_cell_line import read_nucler_total_csv, check_and_remove_duplicated_current_genes, find_overlap

### PARAMETERS
CELL_LINE = "M397"


# read tables
output_file = (project_dir / "02_Clean_data" / "Project_Gene_Annotation_Master.csv")
final_gene_annotation = pd.read_csv(output_file)
clean_dir = project_dir / "02_Clean_data"
master_datasets = load_master_datasets(clean_dir)
print(f"master_datasets: {master_datasets.keys()}")

# read cell line
m397_nuclear_pc, m397_total_pc = read_nucler_total_csv(final_gene_annotation, master_datasets, CELL_LINE)

# check duplicated Current_Gene
collision_gene_qc_filename = project_dir / "05_QC" /"M397_historical_symbol_mapping_collisions.csv"

m397_nuclear_pc_hc, m397_total_pc_hc = check_and_remove_duplicated_current_genes(m397_nuclear_pc, m397_total_pc, collision_gene_qc_filename)
m397_nuclear_pc_hc.to_csv(f"results/m397_nuclear_pc_hc.csv", index=False)
m397_total_pc_hc.to_csv(f"results/m397_total_pc_hc.csv", index=False)


# find genes overlapped between nuclear and total
find_overlap(m397_nuclear_pc_hc, m397_total_pc_hc)
