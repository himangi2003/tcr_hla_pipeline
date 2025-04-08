# test.py
from config import METADATA_FILE, DATA_FILE, DF_HLA_FILE, REP_FOLDER, TRAIN_FOLDER,PTID_SPLIT_FILE
from utils.io_utils import*
from utils.filter import*
from utils.preprocess import*
from utils.train_framework import*
from utils.test_framework import*
import os
import shutil
import pandas as pd
import os
import shutil
import pandas as pd

def run_testing():
    metadata_df = load_metadata(METADATA_FILE)
    data_df = load_data(DATA_FILE)
    DF_HLA = load_hla_data(DF_HLA_FILE)
    
    if os.path.exists(TEST_FOLDER):
        shutil.rmtree(TEST_FOLDER)
    os.makedirs(TEST_FOLDER)

    v = VfamCDR3(project_folder=TEST_FOLDER, input_zfile=REP_FOLDER, cpus=4)
    filenames = v.get_raw_files()

    interaction_df = pd.read_csv("outputs/train_visit_hla_interaction.csv")
    interaction_df = interaction_df[interaction_df['Error'] == False]
    query_df = interaction_df.rename(columns={'TCR': 'vfamcdr3'})

    auc_df = test(DF_HLA, metadata_df, TEST_PTIDS, filenames, query_df, v, edit_type = "edit1")
    auc_df.to_csv("outputs/test_auc.csv", index=False)
