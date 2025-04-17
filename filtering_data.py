from config import*
from utils.io_utils import*
from utils.filter import*
from utils.preprocess import*
import os
import shutil
import pandas as pd
from tqdm import tqdm


def run_filtering(min_num_ptids, only_novel):
    # run filetring from filter.py filtering framework and save outputs 
    
    metadata_df = load_metadata(METADATA_FILE)
    data_df = load_data(DATA_FILE)
    DF_HLA = load_hla_data(DF_HLA_FILE)

    group_1, group_2 = split_ptids(data_df, PTID_SPLIT_FILE)
    train_ptids = group_1['ptid']
    test_ptids = group_2['ptid']



    train_input = process_amino_acid_counts(
        group_1['ptid'].unique(), data_df, REP_FOLDER,
        min_num_ptids = min_num_ptids,
        only_novel = only_novel
    )

    # Save ptid splits to outputs/
    group_1.to_csv("outputs/train_ptids.csv", index=False)
    group_2.to_csv("outputs/test_ptids.csv", index=False)
    # Save list of tcrs for training inputs
    train_input.to_csv("outputs/train_input.csv", index=False)