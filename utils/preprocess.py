#extract sample dataframe
import os
import zipfile
import pandas as pd

def df_extract(sample_file_name, sample_REP_FOLDER_path, sep='\t'):
    """Extracts TCR data reading pre and post samples from a zip archive folder."""
    keep_cols = ['nucleotide', 'aminoAcid', 'count (templates/reads)', "sequenceStatus", 'frequencyCount (%)',
                 'cdr3Length', 'vMaxResolved', 'vGeneName', 'vGeneNameTies', 'jGeneName', 'jGeneNameTies', 'templates_totals']
    
    if not os.path.exists(sample_REP_FOLDER_path):
        return pd.DataFrame()
    
    with zipfile.ZipFile(sample_REP_FOLDER_path) as z3:
        if sample_file_name not in z3.namelist():
            return pd.DataFrame()
        
        with z3.open(sample_file_name) as f:
            df = pd.read_csv(f, sep=sep, usecols=lambda col: col in keep_cols)
            df['vGeneName'] = df['vGeneName'].fillna(df['vGeneNameTies'].str.split(',').str[0])
            df['jGeneName'] = df['jGeneName'].fillna(df['jGeneNameTies'].str.split(',').str[0])
            df = df[df["sequenceStatus"] == "In"]
            df['V'] = df['vGeneName'].str.split("-").str[0].str.replace('TCRB', '', regex=True)
            df['vfamcdr3'] = df['V'] + df['aminoAcid']
            df = df.rename(columns={'count (templates/reads)': 'counts'})
            return df.dropna(subset=['vfamcdr3'])
#sample_mapping
def add_vaccine_sample_ptid_HLA_info_to_df(df, vaccine_meta_DATA_FILE_path):
    m = vaccine_meta_DATA_FILE_path
    m = m.loc[m['sample_selection_cat'] == "B"]
    list_of_ptid = m['pubid'].unique()
    data_dict = {}

    for ptid in list_of_ptid:
        pre_vax_sample = m.loc[(m['pubid'] == ptid) & (m['timepoint'] == 'Pre-vax, pre-inf'), 'sample_name'].iloc[0] \
                            if not m.loc[(m['pubid'] == ptid) & (m['timepoint'] == 'Pre-vax, pre-inf')].empty else None

        post1_vax_sample = m.loc[(m['pubid'] == ptid) & (m['timepoint'] == 'Post-vax, pre-inf'), 'sample_name'].iloc[0] \
                            if not m.loc[(m['pubid'] == ptid) & (m['timepoint'] == 'Post-vax, pre-inf')].empty else None

        if pre_vax_sample:
            pre_vax_sample = pre_vax_sample.replace("re-extraction", "").strip()
        if post1_vax_sample:
            post1_vax_sample = post1_vax_sample.replace("re-extraction", "").strip()

        if pre_vax_sample or post1_vax_sample:
            data_dict[ptid] = {'pre_df': pre_vax_sample, 'post1_df': post1_vax_sample}

    pre_vax_mapping = {samples['pre_df']: ptid + "_pre" for ptid, samples in data_dict.items() if samples['pre_df']}
    post_vax_mapping = {samples['post1_df']: ptid + "_post" for ptid, samples in data_dict.items() if samples['post1_df']}

    df['ptid_pre_df'] = df.index.map(pre_vax_mapping).fillna("")
    df['ptid_post1_df'] = df.index.map(post_vax_mapping).fillna("")
    df['sample_name'] = df.index

    df['ptid_info'] = df['ptid_pre_df'] + df['ptid_post1_df']

    df['ptid'] = df['sample_name'].apply(lambda x: pre_vax_mapping.get(x, "").split("_")[0] if x in pre_vax_mapping else post_vax_mapping.get(x, "").split("_")[0])

    df.drop(columns=['ptid_pre_df', 'ptid_post1_df'], inplace=True)

    return df.fillna(0)

#filter HLA sample for train_ptids
def filter_hla_samples(DF_HLA, m, train_ptids, filenames):
    hla_patient_data_test = add_vaccine_sample_ptid_HLA_info_to_df(DF_HLA, m)
    hla_patient_data_test = hla_patient_data_test[hla_patient_data_test["ptid"].isin(train_ptids)]
    hla_patient_data_test = hla_patient_data_test.sort_values("ptid_info").set_index("ptid_info")
    sample_to_ptid_info_dict = {row["sample_name"]: index for index, row in hla_patient_data_test.iterrows()}
    filtered_filenames = [filename for filename in filenames if any(key in filename for key in sample_to_ptid_info_dict)]
    return filtered_filenames, hla_patient_data_test