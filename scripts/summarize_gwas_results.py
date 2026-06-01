'''
    Script for generating Summary from SnakeGenium Results
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
import os

input_genes = snakemake.input.genes
input_snps = snakemake.input.snps
input_unitigs = snakemake.input.unitigs
input_annotated_unitigs = snakemake.input.annotated_unitigs
input_annotated_snps = snakemake.input.annotated_snps
input_roary = snakemake.input.roary
output_summary = snakemake.output.summary_tsv

def load_and_tag(filepath, variant_type):
    '''
    This function reads a Pyseer output file, extracts the vital stats,
    standardizes the column names, and tags the row with its source type.
    '''
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        try:
            df = pd.read_csv(filepath, sep="\t")
            pval_col = 'adjusted_pvalue' if 'adjusted_pvalue' in df.columns else 'lrt-pvalue'
            cols_to_keep = ['variant', pval_col, 'beta', 'af']
            available_cols = [c for c in cols_to_keep if c in df.columns]
            df_filtered = df[available_cols].copy()
            df_filtered['variant_type'] = variant_type
            df_filtered = df_filtered.rename(columns={pval_col: 'p-value'})
            return df_filtered
        
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
    
    return pd.DataFrame()

def extract_gene_name(attr_str):
    '''
    Extract gene names from the annotations
    In standard GFF files, Column 8 holds attributes (e.g., ID=x;Name=y;...).
    This looks for 'Name=', 'gene=', or 'ID=' and extracts the value.
    '''
    attrs = dict(
        item.split("=", 1)
        for item in str(attr_str).split(";")
        if "=" in item
    )
    return attrs.get("Name") or attrs.get("gene") or attrs.get("ID") or "Unknown_Gene"

df_genes = load_and_tag(input_genes, "Gene")
df_snps = load_and_tag(input_snps, "SNP")
df_unitigs = load_and_tag(input_unitigs, "Unitig")

if not df_unitigs.empty and os.path.exists(input_annotated_unitigs) and os.path.getsize(input_annotated_unitigs) > 0:
    try:
        ann_u_df = pd.read_csv(input_annotated_unitigs, sep="\t", header=None)
        ann_u_df['gene_annotation'] = ann_u_df[8].apply(extract_gene_name)
        df_unitigs['unitig_id'] = [f"unitig_{i}" for i in df_unitigs.index]
        ann_u_df = ann_u_df.merge(df_unitigs[['variant', 'unitig_id']], left_on=12, right_on='unitig_id', how='left')
        ann_mapping = ann_u_df.groupby('variant')['gene_annotation'].apply(lambda x: ', '.join(x.unique())).reset_index()
        df_unitigs = pd.merge(df_unitigs, ann_mapping, on='variant', how='left')
        df_unitigs['gene_annotation'] = df_unitigs['gene_annotation'].fillna("Intergenic/Unmapped")
        df_unitigs = df_unitigs.drop(columns=['unitig_id'])
        
    except Exception as e:
        print(f"Error parsing annotations: {e}")
        df_unitigs['gene_annotation'] = "Annotation Error"
        
elif not df_unitigs.empty:
    df_unitigs['gene_annotation'] = "Unmapped"

if not df_snps.empty and os.path.exists(input_annotated_snps) and os.path.getsize(input_annotated_snps) > 0:
    try:
        ann_s_df = pd.read_csv(input_annotated_snps, sep="\t", header=None)
        ann_s_df['gene_annotation'] = ann_s_df[8].apply(extract_gene_name)
        ann_s_map = ann_s_df.groupby(12)['gene_annotation'].apply(lambda x: ', '.join(x.unique())).reset_index()
        df_snps = pd.merge(df_snps, ann_s_map, left_on='variant', right_on=12, how='left')
        df_snps['gene_annotation'] = df_snps['gene_annotation'].fillna("Intergenic")
        df_snps = df_snps.drop(columns=[12])

    except Exception as e:
        print(f"Error parsing SNP annotations: {e}")
        df_snps['gene_annotation'] = "Annotation Error"

elif not df_snps.empty:
    df_snps['gene_annotation'] = "Intergenic"

if not df_genes.empty:
    if os.path.exists(input_roary):
        try:
            roary_df = pd.read_csv(input_roary, low_memory=False)
            roary_map = roary_df.set_index('Gene')['Annotation'].to_dict()
            df_genes['gene_annotation'] = df_genes['variant'].map(roary_map).fillna(df_genes['variant'])
        except Exception as e:
            print(f"Error reading Roary CSV: {e}")
            df_genes['gene_annotation'] = df_genes['variant']
    else:
        df_genes['gene_annotation'] = df_genes['variant']

valid_dfs = [df for df in [df_genes, df_snps, df_unitigs] if not df.empty]

if valid_dfs:
    combined_df = pd.concat(valid_dfs, ignore_index=True)
else:
    combined_df = pd.DataFrame()

os.makedirs(os.path.dirname(output_summary), exist_ok=True)

if not combined_df.empty:
    cols_order = ['variant_type', 'gene_annotation', 'variant', 'p-value', 'beta', 'af']
    final_cols = [c for c in cols_order if c in combined_df.columns]
    combined_df = combined_df[final_cols].sort_values(by='p-value')
    combined_df.to_csv(output_summary, sep="\t", index=False)
    
else:
    pd.DataFrame(columns=['variant_type', 'gene_annotation', 'variant', 'p-value', 'beta', 'af']).to_csv(output_summary, sep="\t", index=False)