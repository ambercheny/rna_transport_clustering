

# Clustering methods

### 06.py
1. Scenario
    - A gene that changes strongly: [-1.5, +1.5, 0.0, 0.0] (large vector)
    - A gene that changes moderately: [-0.15, +0.15, 0.0, 0.0] (small vector)
2. Using 06.py, clustering would separate these two genes, even though they have the same shape, but they have different magnitude of change

### 07.py
1. This script solves the problem in 06.py by normalizing every gene to be a cector of 1 (unit vector)
2. Then, when clustering, genes that have the same shape will be clustered together