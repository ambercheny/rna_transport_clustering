# Temporal Dynamics of Nuclear vs. Whole-Cell RNA in Targeted Cancer Therapy

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Bioinformatics](https://img.shields.io/badge/Domain-Computational%20Biology-teal.svg)]()

## 📌 Overview

During targeted cancer therapy (e.g., BRAF/MAPK inhibition in melanoma), cancer cells undergo rapid transcriptional adaptation followed by long-term phenotypic rewiring to establish drug tolerance. 

Standard RNA-seq measures the **accumulated whole-cell (total) pool**, which buffers and conceals rapid, transient transcriptional pulses. This repository contains the end-to-end computational pipeline to **decouple acute nuclear transcription from long-term whole-cell RNA fate**:

1. **Nuclear RNA (0h – 72h)**: Single-nucleus / nuclear-fraction RNA capturing immediate, nascent transcriptional surges upon drug perturbation.
2. **Total RNA (Day 0 – Day 38 + Withdrawal)**: Whole-cell RNA-seq tracking long-term drug tolerance, persister state adaptation, and post-withdrawal dynamics.

---

## 🧬 Key Biological Highlights

* **The 12-Hour Emergency Pulse Program**: Unbiased clustering of early nuclear dynamics identifies a unified acute survival program (Cluster 0) co-inducing iron-sequestration genes (`FTH1`, `FTL`), proteotoxic chaperones (`HSPB1`), and lineage regulators (`SOX10`).
* **Cross-Compartment Fate Divergence**:
  * **Transient Emergency Responders** (`HSPB1`, `SERF2`): Sharp 12h nuclear pulse that burns out and returns to baseline in total RNA as acute stress subsides.
  * **Permanent Metabolic Remodelers** (`FTH1`, `FTL`, `FN1`, `ITGB1`): Early nuclear induction translating into permanent, steady-state whole-cell accumulation that forms the drug-tolerant persister phenotype.
* **Lineage De-differentiation**: Long-term tracking reveals an irreversible collapse of the master melanocytic lineage factor **`SOX10`** and **`MITF`** by Day 38, accompanied by upregulation of mesenchymal markers (`FN1`).

---

## 📂 Repository Structure
* Start with 03.py to analyze different cell line

```text
├── 01_Data_preparation.ipynb          # End-to-end interactive exploration, QC & validation
├── 01_clean_data_to_msater_csv.py     # Parses raw Excel workbooks into clean master tables
│   └── utils_01_clean_data.py         # Workbook loading & schema parsing helpers
├── 02_gencode_mapping.py              # HGNC / GENCODE gene symbol standardization
│   └── utils_02_genecode.py           # Collision resolution & alias lookup utilities
├── 03_prepare_cell_line.py            # Extracts & structures matrices per cell line
│   └── utils_03_cell_line.py          # Cell-line filtering & structuring utilities
├── 04_qc_descriptive.py               # Missingness, distribution & correlation QC
├── 05_qc_expression.py                # Expression scale diagnostics & threshold sensitivity
├── 06_clustering.py                   # Baseline K-means clustering on centered trajectories
├── 07_clustering_by_shape.py          # L2 unit-normalized (cosine-equivalent) shape clustering
├── 08_ts_clustering.py                # Time-series clustering & multi-panel trajectory PDFs
│   └── utils_08_ts_cluster.py         # Multi-panel "spaghetti plot" plotting routines
├── 09_multiple_k_means.py             # Multi-k comparison & stability iterations
│   
├── read_clusters.ipynb                # Interactive inspection of cluster outputs
├── unit_normalize_gene_vector.py      # L2-norm spherical normalization module
├── utils_clustering.py                # Core clustering evaluation & matrix preparation
├── results/                           # Output cluster assignments (CSV), centroids, and figures
└── README.md
```

---

## 🔬 Computational Methodology

### 1. Gene Harmonization & Collision Resolution
* Multi-source historical gene symbols standardized to current official **HGNC / GENCODE v50** identities.
* Rare many-to-one alias collisions conservative-filtered to ensure a high-confidence one-to-one gene universe (15,738 nuclear genes; 18,762 total genes; 15,518 shared genes).

### 2. Low-Expression Filtering Sensitivity
* Unsupervised threshold sensitivity analysis testing candidate cutoffs (0.001 to 1.0).
* Operational cutoff selected at **`Max_expression > 0.01`**, safely removing ~25% dead-weight noise while mathematically guaranteeing that no dynamically responsive genes are discarded (`Excluded_DR_max < 0.01`).

### 3. Trajectory Normalization (Shape vs. Amplitude)
* **Log1p Transform**: `log1p(x) = ln(1 + x)` variance stabilization to tame high-expression outliers.
* **Gene-Wise Mean-Centering**: `x - mean(x)` aligns all genes to baseline 0 without distorting true response amplitudes.
* **L2 Unit-Normalization**: Scaling vectors to `||x||₂ = 1`, enabling spherical K-means (cosine similarity) to group genes strictly by temporal shape regardless of abundance.

---

## 📝 Clustering Methods & Practical Notes

### 06.py
1. Scenario
    - A gene that changes strongly: [-1.5, +1.5, 0.0, 0.0] (large vector)
    - A gene that changes moderately: [-0.15, +0.15, 0.0, 0.0] (small vector)
2. Using 06.py, clustering would separate these two genes, even though they have the same shape, but they have different magnitude of change

### 07.py
1. This script solves the problem in 06.py by normalizing every gene to be a cector of 1 (unit vector)
2. Then, when clustering, genes that have the same shape will be clustered together

### 08.py
1. Use `tslearning`
---

### 4. Cluster Validation & Stability
* Quantitative cluster evaluation using **Silhouette**, **Calinski-Harabasz**, **Davies-Bouldin**, and **Inertia**.
* Multi-seed stability validation across 10 random seeds confirming near-deterministic convergence (**Adjusted Rand Index > 0.91**).

---

## 🚀 Quickstart

### Prerequisites

```bash
# Clone the repository
git clone https://github.com/your-username/nuclear-vs-total-rna-dynamics.git
cd nuclear-vs-total-rna-dynamics

# Install dependencies
pip install numpy pandas scikit-learn tslearn scipy matplotlib seaborn openpyxl
```

### Running the Analysis

```bash
# 1. Clean raw Excel data and standardize gene annotations
python 01_clean_data_to_msater_csv.py
python 02_gencode_mapping.py

# 2. Extract cell-line datasets and run QC diagnostics
python 03_prepare_cell_line.py
python 04_qc_descriptive.py
python 05_qc_expression.py

# 3. Trajectory clustering & model evaluation
python 06_clustering.py           # Baseline centered K-means
python 07_clustering_by_shape.py  # L2 unit-normalized shape K-means
python 08_ts_clustering.py        # Time-series clustering & multi-panel PDF generation
```

---

## 📊 Output Deliverables

* **`k_evaluation.csv`**: Goodness-of-fit metric table across `k` in [2, 10].
* **`tslearn_gene_clusters_*_k{n}.csv`**: Gene-level cluster assignments.
* **`cluster_results_*_k{n}.pdf`**: Multi-panel visualization showing all individual gene traces in gray with bold red centroid trajectory curves.

---

