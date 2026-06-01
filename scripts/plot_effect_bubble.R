
# Script for generating Effect Bubble Plot from SnakeGenium Results
# Copyright (C) 2026  Hector H. Furini

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
})

input_summary <- snakemake@input[["summary"]]
output_plot   <- snakemake@output[["plot"]]
antibiotic    <- snakemake@wildcards[["antibiotic"]]

df <- read.delim(input_summary, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

df_unitigs <- df %>% 
    filter(variant_type == "Unitig" & 
           gene_annotation != "Intergenic/Unmapped" & 
           gene_annotation != "Unmapped" & 
           gene_annotation != "Annotation Error")

if(nrow(df_unitigs) > 0) {
    
    agg_df <- df_unitigs %>%
        group_by(gene_annotation) %>%
        summarise(
            max_log_p = max(-log10(`p-value` + 1e-300), na.rm=TRUE),
            avg_effect = mean(abs(beta), na.rm=TRUE),
            avg_maf = mean(ifelse(af > 0.5, 1 - af, af), na.rm=TRUE),
            n_unitigs = n(),
            .groups = 'drop'
        ) %>%
        arrange(desc(avg_effect)) %>%
        mutate(is_top_effect = row_number() <= 15) %>%
        arrange(desc(n_unitigs))

    title_name <- paste(tools::toTitleCase(antibiotic), "Resistance")

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

    ggsave(output_plot, plot = p, width = 9, height = 6, dpi = 300)
    
} else {
    png(output_plot, width = 800, height = 600)
    plot.new()
    text(0.5, 0.5, "No mapped unitig data available to plot")
    dev.off()
}
