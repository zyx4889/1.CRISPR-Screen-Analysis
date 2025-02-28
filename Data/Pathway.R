install.packages("pathfindR")
install.packages("pak")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("ReactomePA")
install.packages("hrbrthemes")


# 1. Load Packages ------------------------------------------------------------
library(pathfindR)       # Pathway enrichment (v1.6.0)
library(dplyr)           # Data manipulation (v1.1.0)
library(ggplot2)         # Visualization (v3.4.0)
library(stringr)         # String operations (v1.5.0)

# 2. Load Input ---------------------------------------------------------------
input_df <- read.csv ("RA_input_synergy.csv")


# 3. Run Multi-Database Analysis ----------------------------------------------
databases <- c("KEGG", "Reactome", "GO-All", "BioCarta")

for(db in databases) {
  output <- run_pathfindR(
    input_df,
    gene_sets = db,
    pin_name_path = "Biogrid",    # Protein interaction network
    p_val_threshold = 0.05        # Adjusted p-value cutoff
  )
  write.csv(output, paste0("", db, "_results.csv"), row.names = FALSE)
}

# 1. Load Results & Preprocess -----------------------------------------------
output_df <- read.csv("GO-All_results.csv") %>% 
  mutate(
    log_p = -log10(lowest_p),
    Gene_Count = str_count(Down_regulated, ',') + 1
  )

# 2. Select Top 15 Pathways --------------------------------------------------
top15 <- output_df %>% 
  arrange(desc(Fold_Enrichment)) %>% 
  head(15)

# 3. Generate Dot Plot -------------------------------------------------------
ggplot(data = top15, aes(x = Fold_Enrichment, 
                         y = reorder(Term_Description, Fold_Enrichment), 
                         color = log_p, 
                         size = freq)) + 
  geom_point() +
  scale_color_gradient2(low = "darkgreen", mid = "seagreen3", high = "indianred2", 
                        midpoint = median(top15$log_p)) +
  scale_size_continuous(range = c(2, 10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "sans", size = 12),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_blank(),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 3  # Adjust this value to change the plot's aspect ratio
  ) +
  labs(
    y = NULL,
    x = "Fold Enrichment",
    color = "-log10(p-value)",
    size = "Number of Genes"
  )
