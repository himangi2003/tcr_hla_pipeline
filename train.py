# train.py
from config import METADATA_FILE, DATA_FILE, DF_HLA_FILE, REP_FOLDER, TRAIN_FOLDER ,PTID_SPLIT_FILE
from utils.io_utils import*
from utils.filter import*
from utils.preprocess import*
from utils.train_framework import*
from utils.test_framework import*
import os
import shutil
import pandas as pd

def run_training():
    metadata_df = load_metadata(METADATA_FILE)
    data_df = load_data(DATA_FILE)
    DF_HLA = load_hla_data(DF_HLA_FILE)


    group_1, group_2 = split_ptids(data_df, PTID_SPLIT_FILE)
    train_ptids = group_1['ptid']
    test_ptids = group_2['ptid']

    processed_group_1 = process_amino_acid_counts(
        group_1['ptid'].unique(), data_df, REP_FOLDER,
        min_num_ptids=2,
        only_novel=True
    )


    # Setup folder
    if os.path.exists(TRAIN_FOLDER):
        shutil.rmtree(TRAIN_FOLDER)
    os.makedirs(TRAIN_FOLDER)

    v = VfamCDR3(project_folder=TRAIN_FOLDER, input_zfile=REP_FOLDER, cpus=4)
    filenames = v.get_raw_files()

    input_tcrs = processed_group_1["vfamcdr3"].unique()

    tcr_presence_absence, tcr_specific_hla, visit, all_fdr_values = analyze_tcr_hla_association(
        DF_HLA, metadata_df, train_ptids, input_tcrs, v, filenames
    )

    hla_patient_data = add_vaccine_sample_ptid_HLA_info_to_df(DF_HLA, metadata_df)
    hla_patient_data = hla_patient_data[hla_patient_data["ptid"].isin(train_ptids)]
    hla_patient_data = hla_patient_data.sort_values("ptid_info").set_index("ptid_info")

    results = {'hla': [], 'visit_hla': [], 'visit_hla_interaction': []}
    for tcr in input_tcrs:
        tcr_result = process_tcr_glm_net(tcr, tcr_presence_absence, visit, tcr_specific_hla, hla_patient_data)
        for key in results:
            results[key].extend(tcr_result[key])

    for key in results:
        df = pd.concat(results[key], ignore_index=True)
        df.to_csv(f"outputs/train_{key}.csv", index=False)
