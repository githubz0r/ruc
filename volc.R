library(tidyverse)
library(magrittr)
library(ggfittext)
library(gt)
library(cowplot)
library(ggVennDiagram)
library(readxl)
library(broom)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Bt.eg.db)
library(gprofiler2)
library(enrichplot)
library(ggVennDiagram)
library(VennDiagram)
library(ggvenn)

data <- read_excel('C:/Users/lassp/Documents/chembio/omics/data (with all grafics).xlsx', sheet='rawdata')
#sample <- data[1, ]
#milk <- data[2,]
#group <-data[3, ]
#data_values <- data[4:nrow(data), ]
data_transposed <- data %>% t() %>% as_tibble() %>% set_colnames(.[1, ]) %>% dplyr::slice(-1)

data_transposed_long <- data_transposed %>% pivot_longer(cols = -c(id, milk, group), names_to = 'ppm', values_to = 'intensity')
data_transposed_long %<>% mutate(milk_simple = ifelse(milk == 'M1', 'human', 'babycow'), intensity = as.numeric(intensity))
results_df <- data_transposed_long %>%
  group_by(ppm) %>%
  summarise(
    mean_human = mean(2^(intensity[milk_simple == 'human'])),
    mean_babycow = mean(2^(intensity[milk_simple == 'babycow'])),
    t_test_results = list(tidy(t.test(intensity ~ milk_simple, data = cur_group()))),
    .groups = 'drop'
  ) %>% 
  unnest(t_test_results) %>% 
  mutate(padjust = p.adjust(p.value, method = "BH"), 
         neglog10padjust = -log10(padjust),
         log2fc = log2(mean_human/mean_babycow),
         significant = ifelse(padjust < 0.05, 'significant', 'not_significant'))




#test <- data_transposed_long %>% filter(ppm == "8.4990000000000006")
#t.test(test %>% filter(milk_simple == 'human') %>% pluck('intensity'), test %>% filter(milk_simple == 'babycow') %>% pluck('intensity'))




# --- Define cutoffs (if not already defined in the dataframe) ---
# LogFC_Cutoff <- 1.0
# FDR_Cutoff <- 0.05 

volcano_plot <- ggplot(
  data = results_df, 
  aes(x = log2fc, y = neglog10padjust, color = significant)
) +
  # 1. Add threshold lines
  # Horizontal line marks the significance (e.g., FDR < 0.05)
  #geom_hline(yintercept = -log10(FDR_Cutoff), linetype = "dashed", color = "gray50") +
  # Vertical lines mark the fold change magnitude (e.g., |Log2FC| > 1.0)
  #geom_vline(xintercept = c(-LogFC_Cutoff, LogFC_Cutoff), linetype = "dashed", color = "gray50") +
  
  # 2. Add the data points
  geom_point(alpha = 0.6, size = 1.5) +
  
  # 3. Apply custom colors to match the 'Expression' column
  #scale_color_manual(values = c("Up-regulated (FDR < 0.05)" = "red", 
                                #"Down-regulated (FDR < 0.05)" = "blue", 
                                #"Not significant" = "gray70")) +
  
  # 4. Customize labels and theme
  labs(
    title = "Volcano Plot: Differential Intensity Analysis (FDR Corrected)",
    #x = expression("Log"[2]*" Fold Change"), # Use expression() for subscript
    #y = expression("-Log"[10]*" Adjusted P-value (FDR)"),
    color = "Significance Status"
  ) +
  theme_minimal() +
  theme(legend.title = element_text(face = "bold"))

print(volcano_plot)



human_vs_formula <- read_tsv('C:/Users/lassp/Documents/chembio/omics/human_vs_formula.tsv')
human_vs_cow <- read_tsv('C:/Users/lassp/Documents/chembio/omics/human_vs_cow.tsv')
cow_vs_formula <- read_tsv('C:/Users/lassp/Documents/chembio/omics/cow_vs_formula.tsv')

human_vs_formula_significant <- human_vs_formula %>% filter(Significant == '+' & !is.na(`Gene names`))
human_vs_cow_significant <- human_vs_cow %>% filter(Significant == '+' & !is.na(`Gene names`))
cow_vs_formula_significant <- cow_vs_formula %>% filter(Significant == '+')

human_vs_formula_up_proteins <- human_vs_formula %>% filter(Significant == '+' & Difference > 0) %>% pluck('Majority protein IDs')
human_vs_formula_down_proteins <- human_vs_formula %>% filter(Significant == '+' & Difference < 0) %>% pluck('Majority protein IDs')

human_vs_cow_up_proteins <- human_vs_cow %>% filter(Significant == '+' & Difference > 0) %>% pluck('Majority protein IDs')
human_vs_cow_down_proteins <- human_vs_cow %>% filter(Significant == '+' & Difference < 0) %>% pluck('Majority protein IDs')

up_prots <- list(human_vs_form_up = human_vs_formula_up_proteins, human_vs_cow_up = human_vs_cow_up_proteins)
down_prots <- list(human_vs_form_down = human_vs_formula_down_proteins, human_vs_cow_down = human_vs_cow_down_proteins)

ggvenn(up_prots)
ggvenn(down_prots)

human_vs_formula_gene_names <- human_vs_formula_significant %>% 
  pluck('Gene names') %>% sapply(function(x){
    if (str_detect(x, ';'))
      {str_split(x, ';')[[1]][1]}
    else x
    })

human_vs_formula_lfc <- setNames(human_vs_formula_significant$Difference, human_vs_formula_gene_names)# %>% sort(decreasing= T)

# human_vs_cow_gene_names <- human_vs_cow %>% 
#   filter(Significant == '+' & !is.na(`Gene names`)) %>% 
#   pluck('Gene names') %>% sapply(function(x){
#     if (str_detect(x, ';'))
#     {str_split(x, ';')[[1]][1]}
#     else x
#   })

enriched_human_vs_formula <- enrichGO(gene = human_vs_formula_gene_names,
         OrgDb = org.Hs.eg.db,
         keyType = "SYMBOL", # maybe change this to gene
         ont = "ALL",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.2)


cnetplot(enriched_human_vs_formula, 
         foldChange = human_vs_formula_lfc, 
         circular   = TRUE, 
         colorEdge  = TRUE)
# enriched_human_vs_cow <- enrichGO(gene = human_vs_cow_gene_names,
#                                       OrgDb = org.Hs.eg.db,
#                                       keyType = "SYMBOL", # maybe change this to gene
#                                       ont = "ALL",
#                                       pAdjustMethod = "BH",
#                                       pvalueCutoff = 0.05,
#                                       qvalueCutoff = 0.2)

# barplot(enriched_human_vs_formula,
#         x = "GeneRatio",
#         color = "p.adjust",
#         title = "Top 10 GO enrichment human vs formula",
#         showCategory = 10,
#         label_format = 80) + facet_grid(.~.sign)
# 
# barplot(enriched_human_vs_cow,
#         x = "GeneRatio",
#         color = "p.adjust",
#         title = "Top 10 GO enrichment human vs cow",
#         showCategory = 10,
#         label_format = 80
# )
# orthologs <- gorth(query = bovine_ids, 
#                    source_organism = "btaurus", 
#                    target_organism = "hsapiens")
# 
# annotations <- bitr("P02702", 
#                     fromType = "UNIPROT", 
#                     toType   = c("SYMBOL", "ENTREZID", "GENENAME"), 
#                     OrgDb    = orgs.Bt.eg.db)


# to do: 

bradford_standard <- function(y){
  x <- (y - 0.0051)/0.0049
  return(x)
}

#absorbances <- c(52.43, 16.10, 6.51, 8.55, 11.81)
#absorbances %>% sapply(bradford_standard)


vec <- c(8.79168, 8.33693, 7.3783, 15.86289, 6.47803, 15.12316994, 14.40885787, 16.42678366, 10.18828, 10.2051, 8.23595, 8.04687, 7.02888, 13.25466, 9.79628, 9.35447)
qqnorm(vec)
qqline(vec)
ggplot(data=NULL, aes(x=vec))+geom_histogram(bins = 7)
