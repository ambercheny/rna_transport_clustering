import numpy as np
import pandas as pd
import pickle
from scipy.spatial.distance import pdist, squareform

def unit_normalize_gene_vector(X):

    X_centered = X.copy()

    # Euclidean norm of each gene trajectory
    row_norms = np.sqrt((X_centered ** 2).sum(axis=1))

    print("Genes with zero trajectory norm:",
        (row_norms == 0).sum())

    # Avoid division by zero
    X_unit = X_centered.div(
        row_norms.replace(0, np.nan),
        axis=0
    )

    print("Unit-normalized matrix shape:", X_unit.shape)
    print("Missing values:", X_unit.isna().sum().sum())
    print(X_unit.head())

    ### check anchor genes
    anchor_genes = ["FTH1", "FTL", "SPRY4", "ETV4"]
    anchor_unit = X_unit.loc[
        X_unit.index.intersection(anchor_genes)
    ]
    print(f"Check gene expression after normaliation of anchor genes")
    print(anchor_unit.round(4))

    print(f"Check pariwise distance")
    anchor_dist = pd.DataFrame(
        squareform(pdist(anchor_unit, metric="euclidean")),
        index=anchor_unit.index,
        columns=anchor_unit.index
    )
    print(anchor_dist.round(4))

    return X_unit
