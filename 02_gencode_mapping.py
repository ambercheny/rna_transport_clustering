# ============================================================
# Parse GENCODE v50 gene and HGNC annotations 
# ============================================================s


import pandas as pd
from pathlib import Path
import openpyxl
import re
import gzip
import re
import pandas as pd
project_dir = Path.cwd().parent
print(f"project_dir: {project_dir}")

from load_master_datasets import load_master_datasets
from utils_02_genecode import read_gencode, merge_project_with_genecode,\
genecode_mapping_status, check_gencode_unmatched_genes,\
read_hgnc, test_hgnc, build_hgnc_symbol_map, label_with_hgnc_rescue,\
validate_genecode_with_hgnc, get_rescued_biotypes_from_hgnc, build_master_gene_annotation,\
qc_duplicated_gene_rows, build_final_gene_annotation, save_final_gene_annotation,\
save_unresolved_genes


qc_dir = project_dir / "05_QC"
gene_summary = pd.read_csv(qc_dir / "Project_unique_gene_symbols.csv")
clean_dir = project_dir / "02_Clean_data"
master_datasets = load_master_datasets(clean_dir)
print(f"master_datasets: {master_datasets.keys()}")
reference_dir = project_dir / "06_Reference"

### GENCODE v50
# read GENCODE v50
gencode_file = reference_dir / "gencode.v50.annotation.gtf.gz"
print("GENCODE file exists:", gencode_file.exists())
gene_annotation = read_gencode(gencode_file)

# merge genecode with our dataset
mapped_genes, unmapped_genes, project_annotation = merge_project_with_genecode(gene_summary, gene_annotation)

# classify genecode mapping status
mapping_summary = genecode_mapping_status(project_annotation)
unmatched_genes = check_gencode_unmatched_genes(mapping_summary, master_datasets)

### use HGNC for unmatched genes
# read hgnc
hgnc_file = reference_dir / "hgnc_complete_set.txt"
print("HGNC file exists:", hgnc_file.exists())
hgnc = read_hgnc(hgnc_file)

# sanity check with AARS
test_hgnc(hgnc, test_symbol = "AARS")

# preprocessing hgnc file, build historial-symnol map
previous_map, alias_map = build_hgnc_symbol_map(hgnc)

# label GENCODE unmatched genes with resecue method
rescue_df = label_with_hgnc_rescue(hgnc, unmatched_genes,previous_map,alias_map)

# validate HGNC with GENCODE
rescued_validated = validate_genecode_with_hgnc(rescue_df, gene_annotation)

# get rescued biotypes from HGNC
rescued_biotypes = get_rescued_biotypes_from_hgnc(rescued_validated)

# create master gene annotation, combining GENCODE and HGNC
project_gene_master = build_master_gene_annotation(gene_summary, gene_annotation, rescued_validated)
qc_duplicated_gene_rows(project_gene_master, gene_summary)

# create final_gene_annotation
final_gene_annotation = build_final_gene_annotation(project_gene_master, gene_summary)

### save files
# save final_gene_annotation
output_file = (
    project_dir
    / "02_Clean_data"
    / "Project_Gene_Annotation_Master.csv"
)
save_final_gene_annotation(final_gene_annotation, output_file)

# save unresolved genes
qc_file = (
    project_dir
    / "05_QC"
    / "Gene_Annotation_Unresolved_Ambiguous.csv"
)
save_unresolved_genes(final_gene_annotation, qc_file)

### summary
print("final_gene_annotation:",final_gene_annotation.shape)
print("M397 nuclear:", master_datasets["M397_nuclear"].shape)
print("M397 total:",master_datasets["M397_total"].shape)