"""Helpers for output creation & export."""
## ============================= IMPORT =============================
import numpy as np
from scipy import stats

## ============================= RESULTS STORAGE =============================
# Mapping from column key to LaTeX header label
COL_LABELS = {
    "sG":                  "$s_G(\\mathbf X)$",
    "subspace invariance": "Subspace invariance",
    "NMI":                 "NMI",
    "AMI":                 "AMI",
    "nb iter":             "Number of iterations",
    "time/iter (avg)":     "Time/iteration",
    "time/iter (ms)":      "Time/iteration (ms)"
}

def init_results_file(file_name, cols, write_type = 'w'):
    col_spec  = "|l|" + "c|" * len(cols)
    headers   = ["Method"] + [COL_LABELS.get(c, c) for c in cols]
    with open(file_name, write_type) as f:
        f.write("\\begin{table}[H]\n")
        f.write("\\centering\n")
        f.write(f"\\begin{{tabular}}{{{col_spec}}}\n")
        f.write("\\hline\n")
        f.write(" & ".join(headers) + " \\\\\n")
        f.write("\\hline\n")


def end_results_file(file_name, example, nb_repet, table_name):
    ex_txt = ""
    if example == "weighted_sbm":
        ex_txt = "synthetic weighted SBM"
    elif example == "caltech":
        ex_txt = "Caltech-7"
    elif example == "synthetic1":
        ex_txt = "synthetic graph 1"
    with open(file_name, 'a') as f:
        f.write("\\hline \n")
        f.write("\\end{tabular}\n")
        f.write("\\caption{Average results with $95\%$ confidence interval for " + ex_txt + " dataset, on " + str(nb_repet) + " runs of the methods.}\n")
        f.write("\\label{tab:" + table_name + "}\n")
        f.write("\\end{table}")
        
def write_results_file(file_name, df, combined):
    for index, row in df.iterrows():
        values = row.values
        if combined:
            formatted = [index] + list(values)
        else:
            formatted = [index] # this is the name of the method
            for i, val in enumerate(values):
                if i == len(values) - 2:  # second to last is integer (nb of iterations)
                    formatted.append(str(int(val)))
                else:  # float
                    formatted.append(f"{val:.3g}")

        output_str = " & ".join(formatted) + " \\\\\n"

        with open(file_name, 'a') as f:
            f.write(output_str)
            
def resultdf_to_latex(df, file_name, example, nb_repet, combined = False):
    init_results_file(file_name)
    # write every line
    write_results_file(file_name, df, combined)
    end_results_file(file_name, example, nb_repet)
    
def resultdf_to_latex_2tables(df, file_name, example, nb_repet, score_keys, extra_keys, combined = False):
    # write first table with scores
    init_results_file(file_name, score_keys)
    write_results_file(file_name, df[score_keys], combined)
    end_results_file(file_name, example, nb_repet, example + "_results")
    # write second table with nb iter + time
    init_results_file(file_name, extra_keys, 'a') # don't overwrite file but add
    write_results_file(file_name, df[extra_keys], combined)
    end_results_file(file_name, example, nb_repet, example + "_iter_time")

def compute_stats(values, confidence=0.95):
    """Return (mean, std, ci_low, ci_high) for a list of values."""
    n = len(values)
    mean = np.mean(values)
    std = np.std(values, ddof=1)
    if n>1:
        se = stats.sem(values)
        h = se*stats.t.ppf((1+confidence)/2, df=n-1)
    else:
        h = np.nan
    return mean, std, mean - h, mean + h

def mean_ci_results_to_file(df_mean, df_std, df_ci_low, df_ci_upp, SCORE_COLS, confidence_label="95% CI", use_ci=True):
    """
    Returns a single DataFrame with entries formatted as 'mean ± value',
    where value is either the CI half-width or the std.
    """
    if use_ci:
        half_width = (df_ci_upp - df_ci_low) / 2
        spread = half_width
        label = confidence_label
    else:
        spread = df_std
        label = "std"

    combined = df_mean.copy().astype(object)
    for col in SCORE_COLS:
        if col in ["sG", "time/iter (avg)"]:
            str_format = "{:.3f}"
        else:
            str_format = "{:.2f}"
        combined[col] = df_mean[col].map(str_format.format) + " \\textpm " + \
                        spread[col].map(str_format.format)
    print(f"Values reported as mean ± {label}\n")
    return combined
