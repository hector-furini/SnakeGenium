'''
    Filtering Significant Hits Script for SnakeGenium
    Copyright (C) 2026  Hector H. Furini

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import pandas as pd
import subprocess
import os
from statsmodels.stats.multitest import multipletests

input_results = snakemake.input.results
input_patterns = snakemake.input.patterns
output_significant = snakemake.output.significant
params_script = snakemake.params.script
wildcards_type = snakemake.wildcards.type

df = pd.read_csv(input_results, sep="\t")
pval_col = "lrt-pvalue" if "lrt-pvalue" in df.columns else "p-value"

if df.empty or pval_col not in df.columns:
    df.head(0).to_csv(output_significant, sep="\t", index=False)
else:
    if wildcards_type == "unitig" and os.path.exists(input_patterns):
        result = subprocess.run(["python", params_script, input_patterns], capture_output=True, text=True)
        p_val_threshold = 0.05 / len(df) if len(df) > 0 else 1e-5
        for line in result.stdout.splitlines():
            if line.startswith("Threshold:"):
                try:
                    val_str = line.split()[1].strip()
                    p_val_threshold = float(val_str)
                except (IndexError, ValueError) as e:
                    print(f"Warning: Could not parse Threshold string. Defaulting to Bonferroni. Error: {e}")

        sig_df = df[df[pval_col] < p_val_threshold].sort_values(by=pval_col)
    
    else:
        reject, pvals_corrected, _, _ = multipletests(df[pval_col], alpha=0.05, method='fdr_bh')
        sig_df = df[reject].copy()
        sig_df['adjusted_pvalue'] = pvals_corrected[reject]
        sig_df = sig_df.sort_values(by=pval_col)
    sig_df.to_csv(output_significant, sep="\t", index=False)