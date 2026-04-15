import pandas as pd
import subprocess
import os
from statsmodels.stats.multitest import multipletests

# Acess SNAKEMAKE files by name
input_results = snakemake.input.results
input_patterns = snakemake.input.patterns
output_significant = snakemake.output.significant
params_script = snakemake.params.script
wildcards_type = snakemake.wildcards.type

df = pd.read_csv(input_results, sep="\t")
# Checks if "lrt-pvalue" exists
pval_col = "lrt-pvalue" if "lrt-pvalue" in df.columns else "p-value"

if df.empty or pval_col not in df.columns:
    # Return empty dataframe with headers if missing p-value column or empty
    df.head(0).to_csv(output_significant, sep="\t", index=False)
else:
    # Calculate dynamic threshold using pyseer script (according to pyseer tutorial)
    # Calculate for unitigs:
    if wildcards_type == "unitig" and os.path.exists(input_patterns):
        # Run pyseer script:
        result = subprocess.run(["python", params_script, input_patterns], capture_output=True, text=True)
        
        # Fail Safe in case the script fails
        p_val_threshold = 0.05 / len(df) if len(df) > 0 else 1e-5

        # Parse stdout string from the script:
        for line in result.stdout.splitlines():
            if line.startswith("Threshold:"):
                try:
                    # Split string and grab value
                    val_str = line.split()[1].strip()
                    p_val_threshold = float(val_str) # Re-write variable
                except (IndexError, ValueError) as e:
                    print(f"Warning: Could not parse Threshold string. Defaulting to Bonferroni. Error: {e}")

        # Filter the final dataframe for significant values using the pattern threshold
        sig_df = df[df[pval_col] < p_val_threshold].sort_values(by=pval_col)
    
    else:
        # For genes and SNPs (where count_patterns.py script isn't used)
        # We will use Benjamini-Hochberg FDR correction via statsmodels
        
        # 'reject' is a boolean array (True if the variant passes the FDR threshold)
        reject, pvals_corrected, _, _ = multipletests(df[pval_col], alpha=0.05, method='fdr_bh')
        
        # Filter the dataframe keeping only the significant hits
        sig_df = df[reject].copy()
        
        # Add the FDR-adjusted p-values as a new column for clarity in the final report
        sig_df['fdr_adjusted_p-value'] = pvals_corrected[reject]
        
        # Sort by the original p-value
        sig_df = sig_df.sort_values(by=pval_col)

    # Export to TSV
    sig_df.to_csv(output_significant, sep="\t", index=False)