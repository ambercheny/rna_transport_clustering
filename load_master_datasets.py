import pandas as pd
def load_master_datasets(clean_dir): 
    """Loads all saved master CSVs back into a dictionary.""" 
    master_datasets = {}
    for file in sorted(clean_dir.glob("*_master.csv")): 
        name = file.name.replace("_master.csv", "")
        master_datasets[name] = pd.read_csv(file)
    return master_datasets

