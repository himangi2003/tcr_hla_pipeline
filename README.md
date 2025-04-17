## Function Reference Table

| Function                                      | Stage         | Description                                                                                 |
|----------------------------------------------|---------------|---------------------------------------------------------------------------------------------|
| `df_extract(...)`                             | Preprocessing | Extracts and cleans TCR sequences from zipped repertoire files, creating `vfamcdr3` identifiers. |
| `add_vaccine_sample_ptid_HLA_info_to_df(...)` | Preprocessing | Maps sample names to PTID and visit labels (`pre`, `post`) in the HLA matrix.              |
| `process_patient(...)`                        | Filtering     | Processes individual patients to identify novel expanding TCRs with count > 1.              |
| `process_amino_acid_counts(...)`              | Filtering     | Aggregates TCRs across patients; retains those seen in ≥ `min_num_ptids`.                   |
| `split_ptids(...)`                            | Filtering     | Splits patients into reproducible training and testing groups using a fixed seed.           |
| `run_filtering(...)`                          | Filtering     | Master function to perform filtering, generate train/test splits, and save TCR input.       |
| `analyze_tcr_hla_association(...)`            | Training      | Builds presence matrix, runs Fisher’s test, and identifies TCR-HLA associations.           |
| `process_tcr_glm_net(...)`                    | Training      | Fits logistic regression models (`~ HLA`, `~ visit + HLA`, `~ visit:HLA`) for each TCR.     |
| `run_training(...)`                           | Training      | Full training pipeline using filtered TCRs, saving GLM results to CSV.                     |
| `compute_auc(...)`                            | Testing       | Computes AUC for each individual HLA subgroup using logistic regression.                    |
| `compute_auc_all_hla(...)`                    | Testing       | Computes AUC across all HLA-positive patients collectively.                                 |
| `test(...)`                                   | Testing       | End-to-end evaluation: prepares test set, builds matrix, computes AUCs.                     |
| `filter_hla_samples(...)`                     | Utility        | Filters filenames and HLA data to match specified PTIDs.                                   |
| `load_metadata(...)`                          | I/O Utility    | Loads and filters metadata file for category 'B'.                                          |
| `load_data(...)`                              | I/O Utility    | Loads main patient-level data file (e.g., PTID, timepoints).                               |
| `load_hla_data(...)`                          | I/O Utility    | Loads the HLA matrix CSV and sets `ptid_info` as index.                                     |
