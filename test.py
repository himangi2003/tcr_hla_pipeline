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

# function to run the testing tramework 
def run_testing(edit_type):
    # Load necessary files
    metadata_df = load_metadata(METADATA_FILE)
    data_df = load_data(DATA_FILE)
    DF_HLA = load_hla_data(DF_HLA_FILE)

    # upload the test PTIDS
    TEST_PTIDS = pd.read_csv("outputs/test_ptids.csv")["ptid"]

    # make test folder for tabify
    if os.path.exists(TEST_FOLDER):
        shutil.rmtree(TEST_FOLDER)
    os.makedirs(TEST_FOLDER)

    v = VfamCDR3(project_folder=TEST_FOLDER, input_zfile=REP_FOLDER, cpus=4)
    filenames = v.get_raw_files()

    # list different types of training models used for analysis and 
    # run testing framework on all of these models
    model_types = ['hla', 'visit_hla', 'visit_hla_interaction']
    all_auc_results = []

    for model_type in model_types:
        # input for testing framework == output of training framework
        interaction_df_path = f"outputs/train_{model_type}_{edit_type}.csv"
        # path for testing framework modelling output to be stored
        output_auc_path = f"outputs/test_auc_{model_type}_{edit_type}.csv"

        if os.path.exists(interaction_df_path):
            interaction_df = pd.read_csv(interaction_df_path)
            interaction_df = interaction_df[interaction_df['Error'] == False]
            # get the query df for tabify 
            query_df = interaction_df.rename(columns={'TCR': 'vfamcdr3'})
            query_df = query_df[query_df["P_value"] <= 0.05]


            # run the test framework from utils.test_framework 
            auc_df = test(DF_HLA, metadata_df, TEST_PTIDS, filenames, query_df, v, edit_type=edit_type)
            auc_df["model_type"] = model_type

            # write the outputs to csv
            auc_df.to_csv(output_auc_path, index=False)  # Save each model's AUC separately
            all_auc_results.append(auc_df)
            print(f" AUC results written to: {output_auc_path}")
        else:
            print(f"[Error] File not found: {interaction_df_path}")


