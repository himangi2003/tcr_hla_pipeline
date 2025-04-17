import pandas as pd
import random
from tqdm import tqdm
from utils.preprocess import*


def process_patient(ptid, data_df, REP_FOLDER):
    """
    Processesing TCR list for a single patient identifying novel TCR expanders and counts>1.

    Args:
        ptid (str): Patient ID.
        data_df (pd.DataFrame): data file that contains 91 PTIDs and sample infomation
        REP_FOLDER (str): Directory path where repertoire files are stored.

    Returns:
        pd.DataFrame: Filtered TCRs with columns ['vfamcdr3', 'counts', 'novel_expander'].
                      Returns an empty DataFrame if no valid data is found.
    """
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
    return count_summary[count_summary['counts'] > 1]

def process_amino_acid_counts(ptids, data_df, REP_FOLDER, min_num_ptids, only_novel):
    """
    Agrregate list of TCRs list scross multiple patients PTIDS identifying novel TCR expanders and counts>1.

    Args:
        ptid (str): Patient ID.
        data_df (pd.DataFrame): data file that contains 91 PTIDs and sample infomation
        REP_FOLDER (str): Directory path where repertoire files are stored.
        min_num_ptids (int): Minimum number of unique PTIDs.
        only_novel (bool): If True, only novel expanders are included.


    Returns:
        pd.DataFrame:TCRs list with columns ['vfamcdr3', 'num_ptids'].
    """
    results = [process_patient(ptid, data_df, REP_FOLDER) for ptid in tqdm(ptids, desc="Processing patients")]
    final_df = pd.concat([r for r in results if not r.empty], ignore_index=True)
    
    if only_novel:
        final_df = final_df[final_df['novel_expander'] == 1]
    
    grouped = final_df.groupby('vfamcdr3').agg(num_ptids=('ptid', 'nunique')).reset_index()
    return grouped[grouped['num_ptids'] >= min_num_ptids]

def process_amino_acid_counts_for_tcrs(tcr_list, data_df, REP_FOLDER, min_ptid_count, min_templates, only_novel):
    """Processes counts for specific TCRs."""
    ptids = data_df['ptid'].unique()
    aa_summary = process_amino_acid_counts(ptids, data_df, REP_FOLDER, min_ptid_count, only_novel)
    return aa_summary[aa_summary['vfamcdr3'].isin(tcr_list)]

def split_ptids(data_df, PTID_SPLIT_FILE, seed=42):
    """
    Splits patient IDs (PTIDs) into two groups for training and testing.

    Args:
        data_df (pd.DataFrame): data file that contains 91 PTIDs and sample infomation.
        PTID_SPLIT_FILE (str): Path to CSV file with list of 91 PTIDs and column 'ptid'.

    Returns:
        tuple: Two DataFrames (group_1, group_2) corresponding to the training and testing splits.
    """
    valid_ptids = pd.read_csv(PTID_SPLIT_FILE)["ptid"].unique()
    random.seed(seed)
    shuffled_ptids = list(valid_ptids)
    random.shuffle(shuffled_ptids)
    
    mid = len(shuffled_ptids) // 2
    group_1_ids = set(shuffled_ptids[:mid])
    group_2_ids = set(shuffled_ptids[mid:])
    
    group_1 = data_df[data_df['ptid'].isin(group_1_ids)]
    group_2 = data_df[data_df['ptid'].isin(group_2_ids)]
    return group_1, group_2


