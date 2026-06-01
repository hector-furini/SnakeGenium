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

input_snps = snakemake.input.snps
output_plot = snakemake.output.plot_file

def extract_variant_info(var_str):
    '''
    Extracts the variant info from the  'variant' column
    in the SNP results dataframe.
    parts[0] = Cromossome id
    parts[1] = Bp position
    '''
    parts = str(var_str).rsplit('_', 3)
    if len(parts) == 4:
        return parts[0], parts[1]
    return "Unknown", np.nan

df = pd.read_csv(input_snps, sep="\t")

pval_col = 'lrt-pvalue' if 'lrt-pvalue' in df.columns else 'p-value'

if not df.empty and pval_col in df.columns:
    extracted_info = df['variant'].apply(extract_variant_info).apply(pd.Series)
    extracted_info.columns = ['#CHR', 'BP']
    
    log10_p = -np.log10(df[pval_col] + 1e-300)
    
    phandango_df = pd.DataFrame({
        '#CHR': extracted_info['#CHR'],
        'SNP': '.',                    
        'BP': extracted_info['BP'],
        'minLOG10(P)': log10_p,
        'log10(p)': log10_p,
        'r^2': 0                    
    })
    
    phandango_df = phandango_df.dropna(subset=['BP'])
    phandango_df['BP'] = phandango_df['BP'].astype(int)
    phandango_df.to_csv(output_plot, sep="\t", index=False)
else:
    with open(output_plot, 'w') as f:
        f.write("#CHR\tSNP\tBP\tminLOG10(P)\tlog10(p)\tr^2\n")