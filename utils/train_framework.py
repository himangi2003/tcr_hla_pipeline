import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import statsmodels.stats.multitest as smm
from statsmodels.stats.multitest import multipletests
from patsy import dmatrices
import scipy.stats as stats
import warnings
from utils.filter import*
from utils.preprocess import*


from tcrtest.ui import VfamCDR3
from tcrtest.tabulate import tabify, tabify1, testify
from tcrtest.classify import HLApredict

# Suppress unnecessary warnings
warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore")


def run_testify(args):
    col, tcr_presence_absence, Y, cutoff = args
    results = testify(Dx=tcr_presence_absence, binary_matrix=Y[[col]].astype("float64"))
    results["fdr_corrected_p"] = smm.multipletests(results["p_exact"], method="fdr_bh")[1]
    results["HLA"] = col
    results = results[results["p_exact"] < cutoff].copy()
    return results

def analyze_tcr_hla_association(
    DF_HLA, metadata_df, train_ptids, input_tcrs, v, filenames, 
    cutoff=0.1, only_novel=True, edit_type="edit0"):
    """
    Analyze TCR-HLA associations with optional novel expander filtering and tabify method selection.

    Parameters:
        DF_HLA (pd.DataFrame): HLA info with ptid_info.
        metadata_df (pd.DataFrame): Metadata including ptids and visit info.
        train_ptids (pd.Series): Patient IDs in training set.
        input_tcrs (list or np.ndarray): TCRs to analyze.
        v (VfamCDR3 object): VfamCDR3 parser instance.
        filenames (list): Raw filenames for parsing.
        cutoff (float): P-value cutoff for inclusion.
        novel_expander (bool): Whether to restrict to post-treatment samples.
        edit_type (str): 'edit0' for tabify, 'edit1' for tabify1.

    Returns:
        tcr_presence_absence (DataFrame), tcr_specific_hla (DataFrame),
        visit (Series), all_fdr_values (DataFrame)
    """
    filtered_filenames, hla_patient_data_train = filter_hla_samples(
        DF_HLA, metadata_df, train_ptids, filenames
    )

    v.parse_adaptive_files(checklist=filtered_filenames)

    filelist = [os.path.join(v.outdir_vfamcdr3, x) for x in os.listdir(v.outdir_vfamcdr3)]
    input_tcrs_df = pd.DataFrame(input_tcrs, columns=["vfamcdr3"])

    # Use tabify or tabify1 depending on edit_type
    if edit_type == "edit1":
        X0 = tabify1(query=input_tcrs_df, filelist=filelist, on='vfamcdr3',
                     get_col='productive_frequency', cpus=4)
    else:
        X0 = tabify(query=input_tcrs_df, filelist=filelist, on='vfamcdr3',
                    get_col='productive_frequency', cpus=4)

    # Binary matrix: TCR presence/absence
    tcr_presence_absence = (X0 > 0).astype(int)
    tcr_presence_absence.index = input_tcrs_df["vfamcdr3"]

    # Build HLA matrix and align columns
    Y = DF_HLA.loc[tcr_presence_absence.columns]
    Y.index = Y["ptid_info"]
    tcr_presence_absence.columns = Y.index

    # Apply novel expander filtering logic
    if only_novel:
        Y_testify = Y[Y.index.str.contains('_post')].iloc[:, :-3]
    else:
        Y_testify = Y.iloc[:, :-3]

    tcr_presence_absence_testify = tcr_presence_absence[Y_testify.index]

    print("Running HLA-TCR association tests...")
    with Pool(cpu_count()) as pool:
        all_results = pool.map(run_testify, [
            (col, tcr_presence_absence_testify, Y_testify, cutoff)
            for col in Y_testify.columns
        ])

    tcr_specific_hla = pd.concat(all_results, ignore_index=True).sort_values("fdr_corrected_p")
    all_fdr_values = tcr_specific_hla.copy()
    tcr_specific_hla.rename(columns={"i": "vfamcdr3"}, inplace=True)

    visit = pd.Series(tcr_presence_absence.columns, index=tcr_presence_absence.columns).str.contains("_post").astype(int)

    return tcr_presence_absence, tcr_specific_hla, visit, all_fdr_values







def process_tcr_data(tcr, tcr_presence_absence, visit, tcr_specific_hla, hla_patient_data):
    """
    Exctract TCR data to fit the model.

    Parameters:
    - tcr (str):  TCR .
    - tcr_presence_absence (pd.DataFrame): DataFrame with TCRs as rows and patients as columns (binary presence/absence).
    - visit (list or np.array): Visit data (should match the number of patients).
    - tcr_specific_hla (pd.DataFrame): Mapping between TCRs and HLA associations.
    - hla_patient_data (pd.DataFrame): HLA information for patients.

    Returns:
    - data (pd.DataFrame or None): DataFrame containing 'presence' and 'visit' and HLA (optional) columns.
    """
    
    if tcr not in tcr_presence_absence.index:
        print(f"Skipping TCR {tcr}: Not found in tcr_presence_absence")
        return None, None

    presence = tcr_presence_absence.loc[tcr].to_numpy()

    if len(presence) != len(visit):
        print(f"Skipping TCR {tcr}: Mismatched presence ({len(presence)}) and visit ({len(visit)}) lengths.")
        return None, None

    hla_indices = tcr_specific_hla.loc[tcr_specific_hla["vfamcdr3"] == tcr, "binary"].tolist()
    valid_hla_indices = [hla for hla in hla_indices if hla in hla_patient_data.columns]

    if not valid_hla_indices:
        hla_df = None
        valid_hla_indices = None
    else:
        patients = tcr_presence_absence.columns
        hla_df = hla_patient_data.loc[patients, valid_hla_indices]

    # Construct DataFrame for presence and visit
    data = pd.DataFrame({'presence': presence, 'visit': visit})

    if hla_df is not None and not hla_df.empty:
            hla_df = hla_df.reset_index(drop=True)
            data = pd.concat([data, hla_df], axis=1)

    return data,hla_df,valid_hla_indices


    

def process_tcr_glm_net(tcr, tcr_presence_absence, visit, tcr_specific_hla, hla_patient_data, penalty='l2', C=1.0, solver='liblinear'):
    """
    Fit logistic regression models using sklearn with different solvers and regularization methods.
    Returns:
    - results (dict): Dictionary with keys ['hla', 'visit_hla', 'visit_hla_interaction'], each containing
                      a list of DataFrames (one per coefficient or error).
    """
    data, hla_df, valid_hla_indices = process_tcr_data(
        tcr, tcr_presence_absence, visit, tcr_specific_hla, hla_patient_data
    )

    if data is None or not valid_hla_indices:
        print(f"[ERROR] TCR: {tcr}, Message: No valid HLA data available.")
        error_row = {
            'TCR': tcr,
            'HLA': None,
            'Model': None,
            'Coefficient': None,
            'Odds_Ratio': None,
            'CI_Lower': None,
            'CI_Upper': None,
            'P_value': None,
            'FDR_pval': None,
            'Error': True,
            'Error_Message': "No valid HLA data available."
        }
        return {
            'hla': [pd.DataFrame([error_row])],
            'visit_hla': [pd.DataFrame([error_row])],
            'visit_hla_interaction': [pd.DataFrame([error_row])]
        }

    results = {
        'hla': [],
        'visit_hla': [],
        'visit_hla_interaction': []
    }

    models = {
        'hla': lambda hla: f'presence ~ {hla}',
        'visit_hla': lambda hla: f'presence ~ visit + {hla}',
        'visit_hla_interaction': lambda hla: f'presence ~ {hla} + visit:{hla}'
    }

    for hla in valid_hla_indices:
        for model_type, formula_func in models.items():
            formula = formula_func(hla)
            try:
                y, X = dmatrices(formula, data, return_type='dataframe')

                if not check_binary_class(y):
                    error_message = "Response variable must be binary (0 and 1)."
                    raise ValueError(error_message)

                model_fit = fit_regularized_logit_sklearn(y, X, penalty, C, solver)
                model_result_df = extract_model_results_sklearn(model_fit, tcr, model_type, hla, X, y)
                model_result_df["Error"] = False
                model_result_df["Error_Message"] = None

            except Exception as e:
                error_message = str(e)
                print(f"[ERROR] TCR: {tcr}, HLA: {hla}, Model: {model_type}, Message: {error_message}")

                model_result_df = pd.DataFrame([{
                    'TCR': tcr,
                    'HLA': hla,
                    'Model': model_type,
                    'Coefficient': None,
                    'Odds_Ratio': None,
                    'CI_Lower': None,
                    'CI_Upper': None,
                    'P_value': None,
                    'FDR_pval': None,
                    'Error': True,
                    'Error_Message': error_message
                }])

            results[model_type].append(model_result_df)

    return results




def fit_regularized_logit_sklearn(y, X, penalty, C, solver):
    model = LogisticRegression(penalty=penalty, C=C, solver=solver, max_iter=1000)
    model.fit(X, y.values.ravel())  
    return model

def extract_model_results_sklearn(model, tcr, model_type, hla, X, y):
    results = []
    coefs = model.coef_.flatten()
    feature_names = X.columns

    eps = np.finfo(float).eps
    pred_probs = model.predict_proba(X)[:, 1]
    W = np.diag(pred_probs * (1 - pred_probs) + eps)

    try:
        cov_matrix = np.linalg.inv(X.T @ W @ X)
        std_errors = np.sqrt(np.diag(cov_matrix))
        z_scores = coefs / std_errors
        p_values = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
    except np.linalg.LinAlgError:
        std_errors = np.full_like(coefs, np.nan)
        p_values = np.full_like(coefs, np.nan)

    # Apply FDR correction across all p-values (just this model's coefficients)
    fdr_values = fdr_correction(p_values)

    for coef_name, coef_value, std_err, p_val, fdr_pval in zip(feature_names, coefs, std_errors, p_values, fdr_values):
        if coef_name != "Intercept":
            results.append({
                'TCR': tcr,
                'HLA': hla,
                'Model': model_type,
                'Coefficient': coef_name,
                'Odds_Ratio': np.exp(coef_value),
                'CI_Lower': np.exp(coef_value - 1.96 * std_err),
                'CI_Upper': np.exp(coef_value + 1.96 * std_err),
                'P_value': p_val,
                'FDR_pval': fdr_pval,
            })

    return pd.DataFrame(results)
    
def fdr_correction(p_values):
    _, corrected_pvals, _, _ = multipletests(p_values, method='fdr_bh')
    return corrected_pvals

def check_binary_class(y):
    unique_classes = np.unique(y)
    return len(unique_classes) > 1