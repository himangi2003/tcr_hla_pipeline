import pandas as pd
import random
from tqdm import tqdm
from utils.preprocess import*

def process_patient(ptid, data_df, REP_FOLDER):
    """Processes a patient and returns TCRs with novel expanders and counts >1."""
    specific_ptid_data = data_df[data_df['ptid'] == ptid]
    if specific_ptid_data.empty:
        return pd.DataFrame()
    
    pre_df = df_extract(specific_ptid_data['pre_vax'].iloc[0], REP_FOLDER)
    post1_df = df_extract(specific_ptid_data['post1_vax'].iloc[0], REP_FOLDER)
    
    if pre_df.empty and post1_df.empty:
        return pd.DataFrame()
    
    novel_expanders = set(post1_df['vfamcdr3']) - set(pre_df['vfamcdr3'])
    count_summary = post1_df.groupby('vfamcdr3')['counts'].sum().reset_index()
    count_summary['ptid'] = ptid
    count_summary['novel_expander'] = count_summary['vfamcdr3'].isin(novel_expanders).astype(int)
    return count_summary[(count_summary['counts'] > 1)]

def process_amino_acid_counts(ptids, data_df, REP_FOLDER, 
                              min_num_ptids: int = 3, only_novel: bool = True):
    """Processes amino acid counts for a list of patients with options to filter novel TCRs and by patient count."""
    results = [process_patient(ptid, data_df, REP_FOLDER) for ptid in tqdm(ptids, desc="Processing patients")]
    final_df = pd.concat([r for r in results if not r.empty], ignore_index=True)
    
    if only_novel:
        final_df = final_df[final_df['novel_expander'] == 1]
    
    grouped = final_df.groupby('vfamcdr3').agg(num_ptids=('ptid', 'nunique')).reset_index()
    return grouped[grouped['num_ptids'] >= min_num_ptids]

def process_amino_acid_counts_for_tcrs(tcr_list, data_df, REP_FOLDER,
                                       min_ptid_count: int = 3, min_templates: int = 1,
                                       only_novel: bool = True):
    """Processes amino acid counts for a given list of TCRs with filtering conditions."""
    ptids = data_df['ptid'].unique()
    amino_acid_summary = process_amino_acid_counts(ptids, data_df, REP_FOLDER, 
                                                   min_num_ptids=min_ptid_count,
                                                   only_novel=only_novel)
    return amino_acid_summary[amino_acid_summary['vfamcdr3'].isin(tcr_list)]


def split_ptids(data_df, PTID_SPLIT_FILE, seed=42):
    """
    Splits PTIDs into two random groups using a reproducible random seed.
    
    Parameters:
        data_df (pd.DataFrame): The main dataframe containing a 'ptid' column.
        PTID_SPLIT_FILE (str): Path to CSV file with a column named 'ptid'.
        seed (int): Random seed for reproducibility (default=42).
        
    Returns:
        group_1 (pd.DataFrame): DataFrame for first random group of ptids.
        group_2 (pd.DataFrame): DataFrame for second random group of ptids.
    """
    valid_ptids = pd.read_csv(PTID_SPLIT_FILE, sep=',')["ptid"].unique()
    random.seed(seed)
    shuffled_ptids = list(valid_ptids)
    random.shuffle(shuffled_ptids)
    
    split_index = len(shuffled_ptids) // 2
    ptid_group_1 = set(shuffled_ptids[:split_index])
    ptid_group_2 = set(shuffled_ptids[split_index:])
    
    group_1 = data_df[data_df['ptid'].isin(ptid_group_1)]
    group_2 = data_df[data_df['ptid'].isin(ptid_group_2)]
    
    return group_1, group_2
