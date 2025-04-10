# train.py
from config import*
from utils.io_utils import*
from utils.filter import*
from utils.preprocess import*
from utils.train_framework import*
from utils.test_framework import*
import os
import shutil
import pandas as pd
from tqdm import tqdm


def run_training(edit_type, min_num_ptids, only_novel):
    metadata_df = load_metadata(METADATA_FILE)
    data_df = load_data(DATA_FILE)
    DF_HLA = load_hla_data(DF_HLA_FILE)

    train_ptids = pd.read_csv("outputs/train_ptids.csv")["ptid"]
    input_tcrs_list_file  = pd.read_csv("outputs/train_input.csv")
    input_tcrs = input_tcrs_list_file["vfamcdr3"].unique()


    # Setup folder
    if os.path.exists(TRAIN_FOLDER):
        shutil.rmtree(TRAIN_FOLDER)
    os.makedirs(TRAIN_FOLDER)

    v = VfamCDR3(project_folder=TRAIN_FOLDER, input_zfile=REP_FOLDER, cpus=4)
    filenames = v.get_raw_files()

    

    tcr_presence_absence, tcr_specific_hla, visit, all_fdr_values = analyze_tcr_hla_association(
        DF_HLA, metadata_df, train_ptids, input_tcrs, v, filenames, 
        cutoff=0.1, only_novel=True, edit_type=edit_type)

    hla_patient_data = add_vaccine_sample_ptid_HLA_info_to_df(DF_HLA, metadata_df)
    hla_patient_data = hla_patient_data[hla_patient_data["ptid"].isin(train_ptids)]
    hla_patient_data = hla_patient_data.sort_values("ptid_info").set_index("ptid_info")

    results = {'hla': [], 'visit_hla': [], 'visit_hla_interaction': []}
    for tcr in tqdm(input_tcrs, desc="Running GLM for each TCR"):
        tcr_result = process_tcr_glm_net(tcr, tcr_presence_absence, visit, tcr_specific_hla, hla_patient_data)
        for key in results:
            results[key].extend(tcr_result[key])


    for key in results:
        df = pd.concat(results[key], ignore_index=True)
        df.to_csv(f"outputs/train_{key}_{edit_type}.csv", index=False)
