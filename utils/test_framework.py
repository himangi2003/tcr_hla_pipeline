import os
import zipfile
import pandas as pd
import numpy as np
import random
import warnings
from utils.filter import*
from utils.preprocess import*
from utils.train_framework import*
from statsmodels.stats.multitest import multipletests

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from tqdm import tqdm

from tcrtest.ui import VfamCDR3
from tcrtest.tabulate import tabify
from tcrtest.classify import HLApredict

# Suppress warnings for cleaner output
warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore")

def compute_auc(df, hla_columns, tcr_columns):
    auc_results = []
    
    for hla in hla_columns:
        DF_HLA = df[df[hla] == 1].copy()
        
        if DF_HLA.empty:
            print(f"Skipping {hla} (no valid samples in this subgroup)")
            continue
        
        DF_HLA["TCR_Sum"] = DF_HLA[tcr_columns].sum(axis=1, skipna=True)
        DF_HLA = DF_HLA.dropna(subset=["TCR_Sum"])
        DF_HLA["Visit"] = DF_HLA.index.to_series().apply(lambda x: 1 if "_post" in x else 0)
        
        if len(DF_HLA["Visit"].unique()) < 2:
            print(f"Skipping {hla} (only one class present in y)")
            continue
        
        X = DF_HLA[["TCR_Sum"]]
        y = DF_HLA["Visit"]
        
        model = LogisticRegression()
        model.fit(X, y)
        
        y_prob = model.predict_proba(X)[:, 1]
        auc = roc_auc_score(y, y_prob)
        
        auc_results.append({"HLA_Subgroup": hla, "AUC": auc})
    
    return pd.DataFrame(auc_results)


def test(DF_HLA, m, test_ptids, filenames, query_df, v):
    if 'TCR' in query_df.columns:
        query_df = query_df.rename(columns={"TCR": "vfamcdr3"})

    filtered_filenames, hla_patient_data_test = filter_hla_samples(DF_HLA, m, test_ptids, filenames)
    v.parse_adaptive_files(checklist=filtered_filenames)
    
    filelist = [os.path.join(v.outdir_vfamcdr3, x) for x in os.listdir(v.outdir_vfamcdr3)]

    X0 = tabify(query=query_df, filelist=filelist, on='vfamcdr3', get_col='productive_frequency', min_value=None, cpus=2)
    X0.index = query_df["vfamcdr3"]
    test_df = X0.T
    sample_map = {row['sample_name']: idx for idx, row in hla_patient_data_test.iterrows()}
    test_df = test_df.rename(index={k: sample_map[k] for k in test_df.index if k in sample_map})
    hla_subgr = hla_patient_data_test.columns.to_list()[:-2]
    df_merged = hla_patient_data_test.merge(test_df, left_index=True, right_index=True, how='inner')
    
    # Derive visit status
    df_merged["visit"] = df_merged.index.to_series().apply(lambda x: 1 if "_post" in x else 0)
    auc_df = compute_auc(df_merged, hla_subgr, query_df["vfamcdr3"])
    return auc_df
