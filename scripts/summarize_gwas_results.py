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

# ---------------------------------------------------------
# GRAB VARIABLES DIRECTLY FROM SNAKEMAKE
# ---------------------------------------------------------
input_genes = snakemake.input.genes
input_snps = snakemake.input.snps
input_unitigs = snakemake.input.unitigs
input_annotated_unitigs = snakemake.input.annotated_unitigs
input_annotated_snps = snakemake.input.annotated_snps
input_roary = snakemake.input.roary
output_summary = snakemake.output.summary_tsv


# This function reads a Pyseer output file, extracts the vital stats,
# standardizes the column names, and tags the row with its source type.
def load_and_tag(filepath, variant_type):
    # Strict check: Only attempt to read if the file exists AND isn't completely empty.
    # Prevents EmptyDataError.
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        try:
            df = pd.read_csv(filepath, sep="\t")
            
            # Pyseer outputs 'lrt-pvalue' for mixed models, but 'p-value' for fixed effects.
            # We detect which one is present so we can standardize it.
            pval_col = 'lrt-pvalue' if 'lrt-pvalue' in df.columns else 'p-value'
            
            # These are the core statistical columns we care about for the final report.
            cols_to_keep = ['variant', pval_col, 'beta', 'af']
            
            # Filter the dataframe to only include these columns (avoids KeyErrors if one is missing).
            available_cols = [c for c in cols_to_keep if c in df.columns]
            df_filtered = df[available_cols].copy()
            
            # Tag the data so we know if a row came from a Gene, SNP, or Unitig.
            df_filtered['variant_type'] = variant_type
            
            # Standardize the p-value column name so it stacks perfectly later.
            df_filtered = df_filtered.rename(columns={pval_col: 'p-value'})
            return df_filtered
        
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
    
    # Return an empty dataframe if the file was empty or missing.
    return pd.DataFrame()

# In standard GFF files, Column 8 holds attributes (e.g., ID=x;Name=y;...).
# This looks for 'Name=', 'gene=', or 'ID=' and extracts the value.
def extract_gene_name(attr_str):
    attrs = dict(
        item.split("=", 1) # Splits each item by "="
        for item in str(attr_str).split(";") # For each item in column splited by ";"
        if "=" in item # If "=" in column
    )
    return attrs.get("Name") or attrs.get("gene") or attrs.get("ID") or "Unknown_Gene" # Return

# LOAD SIGNIFICANT HITS
df_genes = load_and_tag(input_genes, "Gene")
df_snps = load_and_tag(input_snps, "SNP")
df_unitigs = load_and_tag(input_unitigs, "Unitig")

# ==========================================
# 1. ANNOTATE UNITIGS
# ==========================================
# PARSE BEDTOOLS GFF ANNOTATIONS & MERGE WITH UNITIGS
# Proceed only if we actually have significant unitigs AND an annotation file.
if not df_unitigs.empty and os.path.exists(input_annotated_unitigs) and os.path.getsize(input_annotated_unitigs) > 0:
    try:
        # Bedtools outputs do not have headers. We load it raw.
        ann_u_df = pd.read_csv(input_annotated_unitigs, sep="\t", header=None)

        # Create a clean 'gene_annotation' column from the messy GFF attributes.
        ann_u_df['gene_annotation'] = ann_u_df[8].apply(extract_gene_name)

        # Create the matching unitig_id in the original stats dataframe to serve as a key
        df_unitigs['unitig_id'] = [f"unitig_{i}" for i in df_unitigs.index]
        
        # Column 12 contains 'unitig_X' ID. Map it to the original sequence.
        ann_u_df = ann_u_df.merge(df_unitigs[['variant', 'unitig_id']], left_on=12, right_on='unitig_id', how='left')
        
        # Now that 'variant' holds the real DNA sequence, group by it
        ann_mapping = ann_u_df.groupby('variant')['gene_annotation'].apply(lambda x: ', '.join(x.unique())).reset_index()

        # Map these cleaned annotations back to our Pyseer stats dataframe.
        df_unitigs = pd.merge(df_unitigs, ann_mapping, on='variant', how='left')
        
        # If a unitig mapped to the genome but didn't hit a gene feature, it's likely intergenic.
        df_unitigs['gene_annotation'] = df_unitigs['gene_annotation'].fillna("Intergenic/Unmapped")
        df_unitigs = df_unitigs.drop(columns=['unitig_id'])
        
    except Exception as e:
        print(f"Error parsing annotations: {e}")
        df_unitigs['gene_annotation'] = "Annotation Error"
        
# If we have unitigs but the bedtools file was empty, mark them as unmapped.
elif not df_unitigs.empty:
    df_unitigs['gene_annotation'] = "Unmapped"
# ==========================================
# 2. ANNOTATE SNPs
# ==========================================
# Check for empty or non-existing files
if not df_snps.empty and os.path.exists(input_annotated_snps) and os.path.getsize(input_annotated_snps) > 0:
    try:
        # Bedtools outputs do not have headers. We load it raw.
        ann_s_df = pd.read_csv(input_annotated_snps, sep="\t", header=None)

        # Create a clean 'gene_annotation' column from the messy GFF attributes.
        ann_s_df['gene_annotation'] = ann_s_df[8].apply(extract_gene_name)
        
        # Map the variant string (Column 12) to the extracted annotation
        ann_s_map = ann_s_df.groupby(12)['gene_annotation'].apply(lambda x: ', '.join(x.unique())).reset_index()
        
        # Map these cleaned annotations back to our Pyseer stats dataframe.
        df_snps = pd.merge(df_snps, ann_s_map, left_on='variant', right_on=12, how='left')
        
        # If a SNP mapped to the genome but didn't hit a gene feature, it's likely intergenic.
        df_snps['gene_annotation'] = df_snps['gene_annotation'].fillna("Intergenic")
        df_snps = df_snps.drop(columns=[12]) # Remove original column

    except Exception as e:
        print(f"Error parsing SNP annotations: {e}")
        df_snps['gene_annotation'] = "Annotation Error"

elif not df_snps.empty:
    df_snps['gene_annotation'] = "Intergenic"

# ==========================================
# 3. ANNOTATE GENES (VIA ROARY)
# ==========================================

# Check for empty DataFrame
if not df_genes.empty:
    if os.path.exists(input_roary): # Check for Roary output
        try:
            roary_df = pd.read_csv(input_roary, low_memory=False)
            roary_map = roary_df.set_index('Gene')['Annotation'].to_dict()
            df_genes['gene_annotation'] = df_genes['variant'].map(roary_map).fillna(df_genes['variant'])
        except Exception as e:
            print(f"Error reading Roary CSV: {e}")
            df_genes['gene_annotation'] = df_genes['variant']
    else:
        df_genes['gene_annotation'] = df_genes['variant']

# ==========================================
# COMBINE AND EXPORT
# ==========================================
# Filter out empty DataFrames before passing to pd.concat
valid_dfs = [df for df in [df_genes, df_snps, df_unitigs] if not df.empty]

if valid_dfs:
    # COMBINE AND EXPORT
    # Stack the three dataframes vertically.
    combined_df = pd.concat(valid_dfs, ignore_index=True)
else:
    combined_df = pd.DataFrame() # Create an empty DataFrame if all inputs were empty

# Ensure the target directory actually exists before saving.
os.makedirs(os.path.dirname(output_summary), exist_ok=True)

if not combined_df.empty:
    # Define a clean, readable column order for the final report.
    cols_order = ['variant_type', 'gene_annotation', 'variant', 'p-value', 'beta', 'af']
    final_cols = [c for c in cols_order if c in combined_df.columns]
    
    # Sort the entire merged dataset by p-value (most significant hits at the top) and save.
    combined_df = combined_df[final_cols].sort_values(by='p-value')
    combined_df.to_csv(output_summary, sep="\t", index=False)
    
else:
    # If there are zero significant hits across all categories, generate an empty file 
    # with the correct headers so downstream processes don't break.
    pd.DataFrame(columns=['variant_type', 'gene_annotation', 'variant', 'p-value', 'beta', 'af']).to_csv(output_summary, sep="\t", index=False)