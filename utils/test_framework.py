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

# Suppress warnings 
warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore")

# function to compute AUC for each HLA subgroup at a time 
def compute_auc(df, hla_columns, tcr_columns):
    auc_results = []
  
    # for modeling visit ~ sum(tcr_presence_absence) for each HLA subgroup 
    # if a valid test ptid samples is associated with this subgroup

    for hla in hla_columns:
        DF_HLA = df[df[hla] == 1].copy()

        if DF_HLA.empty:
            print(f"Skipping {hla} (no valid samples in this subgroup)")
            continue

        # find the X variable from df
        DF_HLA["TCR_Sum"] = DF_HLA[tcr_columns].sum(axis=1, skipna=True)
        DF_HLA = DF_HLA.dropna(subset=["TCR_Sum"])

        # find Y variable from df
        DF_HLA["Visit"] = DF_HLA.index.to_series().apply(lambda x: 1 if "_post" in x else 0)

        # check for single class error
        if len(DF_HLA["Visit"].unique()) < 2:
            print(f"Skipping {hla} (only one class present in y)")
            continue

        X = DF_HLA[["TCR_Sum"]]
        y = DF_HLA["Visit"]

        # train logit model 
        model = LogisticRegression(penalty='l2', C=1.0, solver='liblinear')
        model.fit(X, y)

        # find AUC
        y_prob = model.predict_proba(X)[:, 1]
        auc = roc_auc_score(y, y_prob)

        # save the result
        auc_results.append({
            "HLA_Subgroup": hla,
            "AUC": auc,
            "AUC_Source": "per_HLA"
        })

    return pd.DataFrame(auc_results)


def compute_auc_all_hla(df, hla_columns, tcr_columns):

    # for modeling visit ~ sum(tcr_presence_absence) for all HLA at
    # if a valid test ptid samples is associated with this subgroup

    hla_mask = df[hla_columns].sum(axis=1) > 0
    DF_HLA = df[hla_mask].copy()

    if DF_HLA.empty:
        print("Skipping all HLA (no valid samples with any HLA present)")
        return pd.DataFrame([{
            "HLA_Subgroup": "all_HLA",
            "AUC": None,
            "AUC_Source": "all_HLA_combined"
        }])

    DF_HLA["TCR_Sum"] = DF_HLA[tcr_columns].sum(axis=1, skipna=True)
    DF_HLA = DF_HLA.dropna(subset=["TCR_Sum"])
    DF_HLA["Visit"] = DF_HLA.index.to_series().apply(lambda x: 1 if "_post" in x else 0)

    if len(DF_HLA["Visit"].unique()) < 2:
        print("Skipping all HLA (only one class present in Visit column)")
        return pd.DataFrame([{
            "HLA_Subgroup": "all_HLA",
            "AUC": None,
            "AUC_Source": "all_HLA_combined"
        }])

    X = DF_HLA[["TCR_Sum"]]
    y = DF_HLA["Visit"]

    # train logit model 
    model = LogisticRegression()
    model.fit(X, y)

    # find AUC
    y_prob = model.predict_proba(X)[:, 1]
    auc = roc_auc_score(y, y_prob)

    # save the result
    return pd.DataFrame([{
        "HLA_Subgroup": "all_HLA",
        "AUC": auc,
        "AUC_Source": "all_HLA_combined"
    }])



def test(DF_HLA, m, test_ptids, filenames, query_df, v, edit_type="edit0"):
    """
    Evaluate TCR-HLA association models by computing AUC scores on test data.

    This function processes the HLA typing matrix for a set of test patients (PTIDs), 
    extracts TCR presence/absence from repertoire files, and integrates this with HLA typing 
    and visit information (pre- vs. post-treatment). It then computes AUC scores to evaluate 
    the predictive power of TCRs across:
    - individual HLA subgroups (`compute_auc`)
    - all HLA alleles combined (`compute_auc_all_hla`)

    Parameters:
    ----------
    DF_HLA : pd.DataFrame
        HLA typing matrix with binary values (1/0) indicating presence/absence of specific HLA alleles per sample.
    
    m : pd.DataFrame
        Patient metadata including PTID, sample names, and visit information (pre- and post-treatment).
    
    test_ptids : list or pd.Series
        List of patient IDs (PTIDs) to include in the test set.
    
    filenames : list
        List of raw TCR repertoire file paths for Test PTIDs set.
    
    query_df : pd.DataFrame
        DataFrame of query TCRs with a 'TCR' or 'vfamcdr3' column, derived from model training output.
    
    v : VfamCDR3
        obtained using tabify function.
    
    edit_type : str,("edit0", "edit1")
        Specifies which tabify function to use for generating the TCR matrix:
        - "edit0": uses `tabify` (default behavior)
        - "edit1": uses `tabify1` for alternative processing

    Returns:
    -------
    pd.DataFrame
        Combined AUC scores for:
        - each HLA allele separately (output of `compute_auc`)
        - all HLAs together (output of `compute_auc_all_hla`)
    """
  
   #check query df file obtained from train framework output
   # get list of TCRs in query_df file to use tabify

    if 'TCR' in query_df.columns:
        query_df = query_df.rename(columns={"TCR": "vfamcdr3"})

    # preprocesing and finding repetoire filtered filenames list based on test ptids sample to use tabify function 
    # and  HLA typing matrix for filtered filenames list
    filtered_filenames, hla_patient_data_test = filter_hla_samples(DF_HLA, m, test_ptids, filenames)
    v.parse_adaptive_files(checklist=filtered_filenames)

    filelist = [os.path.join(v.outdir_vfamcdr3, x) for x in os.listdir(v.outdir_vfamcdr3)]

    # use tabify/tabify1
    if edit_type == "edit1":
        X0 = tabify1(query=query_df, filelist=filelist, on='vfamcdr3',
                     get_col='productive_frequency', cpus=4)
    elif edit_type == "edit0":
        X0 = tabify(query=query_df, filelist=filelist, on='vfamcdr3',
                    get_col='productive_frequency', cpus=4)
    else:
        print("Please choose a valid edit type: 'edit0' or 'edit1'")
        X0 = None

    # process tabify output
    X0.index = query_df["vfamcdr3"]
    test_df = X0.T

    # map sample names from filtered file list and test HLA typing matrix and create a test dataframe with a common index 
    sample_map = {row['sample_name']: idx for idx, row in hla_patient_data_test.iterrows()}
    test_df = test_df.rename(index={k: sample_map[k] for k in test_df.index if k in sample_map})

    # find HLA subgroups to be tested 
    hla_subgr = hla_patient_data_test.columns.to_list()[:-2]
    df_merged = hla_patient_data_test.merge(test_df, left_index=True, right_index=True, how='inner')

    # find visit matrix (post = 1 , pre = 0 ) which will be used as dependent variable in the compute AUC models
    df_merged["visit"] = df_merged.index.to_series().apply(lambda x: 1 if "_post" in x else 0)

    # Run both AUC computations
    auc_df_per_hla = compute_auc(df_merged, hla_subgr, query_df["vfamcdr3"])
    auc_df_all_hla = compute_auc_all_hla(df_merged, hla_subgr, query_df["vfamcdr3"])

    return pd.concat([auc_df_per_hla, auc_df_all_hla], ignore_index=True)

