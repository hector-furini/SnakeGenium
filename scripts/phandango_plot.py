'''
    Script for generating Manhattan Plot from SnakeGenium Results
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
import numpy as np

# Import SNAKEMAKE variables
input_snps = snakemake.input.snps
output_plot = snakemake.output.plot_file

df = pd.read_csv(input_snps, sep="\t")

# Determine the correct p-value column based on the model (LMM vs Fixed)
pval_col = 'lrt-pvalue' if 'lrt-pvalue' in df.columns else 'p-value'

if not df.empty and pval_col in df.columns:
    
    # Pyseer SNP variants format: {CHROM}_{POS}_{REF}_{ALT}
    # We use rsplit('_', 3) to safely separate it from the right, 
    # ensuring CHROM names with underscores remain intact.
    def extract_variant_info(var_str):
        parts = str(var_str).rsplit('_', 3)
        if len(parts) == 4:
            return parts[0], parts[1] # Returns (CHROM, POS)
        return "Unknown", np.nan
    
    # Apply the function and split results into two columns
    extracted_info = df['variant'].apply(extract_variant_info).apply(pd.Series)
    extracted_info.columns = ['#CHR', 'BP']
    
    # Calculate -log10(p-value), adding tiny offset to prevent log10(0) infinity
    log10_p = -np.log10(df[pval_col] + 1e-300)
    
    # Build the exact dataframe structure per the Pyseer bash tutorial
    phandango_df = pd.DataFrame({
        '#CHR': extracted_info['#CHR'],  # Now perfectly matches the GFF sequence ID
        'SNP': '.',                      # Dot used as placeholder
        'BP': extracted_info['BP'],
        'minLOG10(P)': log10_p,
        'log10(p)': log10_p,
        'r^2': 0                         # Hardcoded to 0 as required
    })
    
    # Drop unparseable variants and format position as integer
    phandango_df = phandango_df.dropna(subset=['BP'])
    phandango_df['BP'] = phandango_df['BP'].astype(int)
    
    # Phandango requires exact header casing and tab separation
    phandango_df.to_csv(output_plot, sep="\t", index=False)
else:
    # Fallback for empty results
    with open(output_plot, 'w') as f:
        f.write("#CHR\tSNP\tBP\tminLOG10(P)\tlog10(p)\tr^2\n")