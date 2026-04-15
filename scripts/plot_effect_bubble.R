# Load required libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
})

# ---------------------------------------------------------
# GRAB VARIABLES DIRECTLY FROM SNAKEMAKE
# ---------------------------------------------------------
input_summary <- snakemake@input[["summary"]]
output_plot   <- snakemake@output[["plot"]]
antibiotic    <- snakemake@wildcards[["antibiotic"]]

# Read the summary file safely, preventing R from changing "p-value" to "p.value"
df <- read.delim(input_summary, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

# 1. Filter for Unitigs mapped to a known gene
df_unitigs <- df %>% 
    filter(variant_type == "Unitig" & 
           gene_annotation != "Intergenic/Unmapped" & 
           gene_annotation != "Unmapped" & 
           gene_annotation != "Annotation Error")

# Check if we have data to plot
if(nrow(df_unitigs) > 0) {
    
    # 2. Aggregate data by gene and sort
    agg_df <- df_unitigs %>%
        group_by(gene_annotation) %>%
        summarise(
            max_log_p = max(-log10(`p-value` + 1e-300), na.rm=TRUE),
            avg_effect = mean(abs(beta), na.rm=TRUE),
            # Calculate Minor Allele Frequency (MAF)
            avg_maf = mean(ifelse(af > 0.5, 1 - af, af), na.rm=TRUE),
            n_unitigs = n(),
            .groups = 'drop'
        ) %>%
        # Plot sorted by effect size
        arrange(desc(avg_effect)) %>%
        mutate(is_top_effect = row_number() <= 15) %>%
        arrange(desc(n_unitigs))

    # Capitalize the antibiotic name for the title
    title_name <- paste(tools::toTitleCase(antibiotic), "Resistance")

    # 3. Generate the ggplot bubble chart
    p <- ggplot(agg_df, aes(x = avg_effect, y = max_log_p)) +
        geom_point(aes(size = n_unitigs, color = avg_maf), alpha = 0.6) +
        geom_text_repel(aes(label = gene_annotation,
                            fontface = ifelse(is_top_effect, "bold", "plain")), 
                        size = 3,
                        max.overlaps = 20,
                        min.segment.length = 0,
                        box.padding = 0.5,
                        show.legend = FALSE) +

        scale_color_gradient(low = "#132B43", high = "#56B1F7", name = "Average MAF") +
        scale_size_continuous(name = "Number of Unitigs") +
        labs(title = title_name,
             x = "Average Effect Size (|beta|)", 
             y = "Maximum -log10(p-value)") +
        theme_bw() +
        theme(panel.grid.minor = element_blank())

    # Save the plot using the Snakemake output path
    ggsave(output_plot, plot = p, width = 9, height = 6, dpi = 300)
    
} else {
    # Safe fallback if no mapped unitigs exist for this antibiotic
    png(output_plot, width = 800, height = 600)
    plot.new()
    text(0.5, 0.5, "No mapped unitig data available to plot")
    dev.off()
}
