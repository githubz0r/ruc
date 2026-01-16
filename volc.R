library(tidyverse)
library(magrittr)
library(ggfittext)
library(gt)
library(cowplot)
library(ggVennDiagram)
library(readxl)
library(broom)
library(ggplot2)

norm_data <- read_excel('/Users/b246357/Documents/uni/data.xlsx', sheet='normdata')
raw_data <- read_excel('/Users/b246357/Documents/uni/data.xlsx', sheet='rawdata')

#sample <- data[1, ]
#milk <- data[2,]
#group <-data[3, ]
#data_values <- data[4:nrow(data), ]
data_transposed <- norm_data %>% t() %>% as_tibble() %>% set_colnames(.[1, ]) %>% slice(-1)

data_transposed_long <- data_transposed %>% pivot_longer(cols = -c(id, milk, group), names_to = 'ppm', values_to = 'intensity')
#data_transposed_long %<>% filter(intensity > 0) # filter out noise, make a function to make sure every class has at least 3 values
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
         #log2fc = log2(mean_babycow/mean_human),
         significant = case_when(
           padjust < 0.01 & log2fc > 0 ~ "higher abundance in human",
           padjust < 0.01 & log2fc < 0 ~ "lower abundance in human",
           TRUE ~ "not significant"
         ))




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
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
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


#bmrb_m2m3 <- read_csv('/Users/b246357/Documents/uni/data_m2_m3_metabolites.csv')
#bmrb_m2m3 %<>% mutate(
  # The regex pattern:
  # (?<=H:) looks behind for "H:" (but doesn't include it in the match)
  # [0-9.]+ matches any digits or decimal points that follow
#  H_Value_Char = str_extract(`Matching shifts...6`, "(?<=H:)[0-9\\.]+")
#)
