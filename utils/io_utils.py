# tcr_hla_pipeline/utils/io_utils.py
import pandas as pd
import os

def load_metadata(metadata_file: str) -> pd.DataFrame:
    """Loads metadata file and filters for category 'B'."""
    df = pd.read_csv(metadata_file, dtype={'sample_selection_cat': str})
    return df[df['sample_selection_cat'] == 'B'] if not df.empty else pd.DataFrame()

def load_data(data_file: str) -> pd.DataFrame:
    """Loads patient-level data from the given CSV file."""
    return pd.read_csv(data_file, dtype={'ptid': str}) if os.path.exists(data_file) else pd.DataFrame()

def load_hla_data(hla_file: str) -> pd.DataFrame:
    """Loads the HLA matrix with ptid_info index."""
    return pd.read_csv(hla_file, index_col=0) if os.path.exists(hla_file) else pd.DataFrame()
