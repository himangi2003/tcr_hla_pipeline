# test.py
from config import*
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

def run_testing(edit_type):
    metadata_df = load_metadata(METADATA_FILE)
    data_df = load_data(DATA_FILE)
    DF_HLA = load_hla_data(DF_HLA_FILE)
    TEST_PTIDS = pd.read_csv("outputs/test_ptids.csv")["ptid"]

    if os.path.exists(TEST_FOLDER):
        shutil.rmtree(TEST_FOLDER)
    os.makedirs(TEST_FOLDER)

    v = VfamCDR3(project_folder=TEST_FOLDER, input_zfile=REP_FOLDER, cpus=4)
    filenames = v.get_raw_files()

    model_types = ['hla', 'visit_hla', 'visit_hla_interaction']
    all_auc_results = []

    for model_type in model_types:
        interaction_df_path = f"outputs/train_{model_type}_{edit_type}.csv"
        output_auc_path = f"outputs/test_auc_{model_type}_{edit_type}.csv"

        if os.path.exists(interaction_df_path):
            interaction_df = pd.read_csv(interaction_df_path)
            interaction_df = interaction_df[interaction_df['Error'] == False]
            query_df = interaction_df.rename(columns={'TCR': 'vfamcdr3'})

            auc_df = test(DF_HLA, metadata_df, TEST_PTIDS, filenames, query_df, v, edit_type=edit_type)
            auc_df["model_type"] = model_type

            auc_df.to_csv(output_auc_path, index=False)  # Save each model's AUC separately
            all_auc_results.append(auc_df)
            print(f" AUC results written to: {output_auc_path}")
        else:
            print(f"[Error] File not found: {interaction_df_path}")


