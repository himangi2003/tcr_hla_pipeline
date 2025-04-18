## Project Directory Structure

- `tcr-hla-pipeline/`
  - `utils/`  
    - `io_utils.py` – Loads metadata, patient data, and HLA data  
    - `preprocess.py` – preprocesses Rep files and prepares input dataframe and maps sample metadata  
    - `train_framework.py` – model training and result extraction  
    - `filter.py` – Novel TCR detection and patient-level filtering  
    - `test_framework.py` – AUC computation and evaluation logic  
  - `config.py` – Configuration file with paths  
  - `filtering.py` – Runs the TCR filtering step  
  - `train.py` – Runs the training step using filtered TCRs  
  - `test.py` – Evaluates models on test set using AUC metrics  
  - `main.py` – Runs the full pipeline (filter → train → test)  
  - `scripts/`  
    - `run_pipeline.slurm` – SLURM script to run the full pipeline on HPC





## Function Reference Table

| Function                                      | Stage         | Description                                                                                 |
|----------------------------------------------|---------------|---------------------------------------------------------------------------------------------|
| `df_extract(...)`                             | Preprocessing | Extracts and cleans Vfamcdr3 dataframe from zipped repertoire files. 			     |
| `add_vaccine_sample_ptid_HLA_info_to_df(...)` | Preprocessing | Maps sample names to PTID and visit labels (`pre`, `post`) in the HLA typing matrix.              |
| `process_patient(...)`                        | Filtering     | Processes individual patients to identify novel expanding TCRs with count > 1.              |
| `process_amino_acid_counts(...)`              | Filtering     | Aggregates TCRs across patients; retains those seen in ≥ `min_num_ptids`.                   |
| `split_ptids(...)`                            | Filtering     | Splits patients into reproducible training and testing groups using a fixed seed.           |
| `run_filtering(...)`                          | Filtering     | Master function to perform filtering, generate train/test splits, and save TCR input.       |
| `analyze_tcr_hla_association(...)`            | Training      | Builds TCR presence/ absence matrix, runs testify and identifies TCR-HLA associations.           |
| `process_tcr_glm_net(...)`                    | Training      | Fits logistic regression models (`~ HLA`, `~ visit + HLA`, `~ visit:HLA`) for each TCR.     |
| `run_training(...)`                           | Training      | Full training pipeline using filtered TCRs, saving GLM results to CSV.                     |
| `compute_auc(...)`                            | Testing       | Computes AUC for each individual HLA subgroup using logistic regression.                    |
| `compute_auc_all_hla(...)`                    | Testing       | Computes AUC across all HLA-positive patients collectively.                                 |
| `test(...)`                                   | Testing       | testing framework: prepares test set, builds matrix, computes AUCs.                     |
| `run_testing(...)`                            | Testing       | Master function to perform Testing, generate HLA subgroups with AUC information                    |
| `filter_hla_samples(...)`                     | Utility        | Filters filenames and HLA data to match specified PTIDs.                                   |
| `load_metadata(...)`                          | I/O Utility    | Loads and filters metadata file for category 'B'.                                          |
| `load_data(...)`                              | I/O Utility    | Loads main patient-level data file (e.g., PTID, timepoints).                               |
| `load_hla_data(...)`                          | I/O Utility    | Loads the HLA matrix CSV and sets `ptid_info` as index.                                     |
