library(tidyverse)
library(magrittr)

exp_cell_path <- "/Users/b246357/Documents/uni/expcell"
w1_cell_area <- read_csv(file.path(exp_cell_path, "exp1_w1_cell_areas.csv"))
w1_cell_area %<>% mutate(cell_area = `Cell Area`/ncell) # average each mass of cells
inch_to_um <- (300^2)/(2.96^2) # goes from inches^2 (imagej default) to µM^2
w1_cell_area %<>% mutate(cell_area = inch_to_um * cell_area)
w1_cell_area %>% ggplot(aes(x= mM_NaBu, y = cell_area, col = cell_line))+geom_jitter()

w1_cell_area_summary_line <- w1_cell_area %>% group_by(cell_line, Well, mM_NaBu) %>% 
  summarise(mean_area = mean(cell_area),.groups = "drop")

w1_cell_area_summary_line <- w1_cell_area_summary_line %>%
  mutate(
    well_letter = str_extract(Well, "^[A-Z]"),
    well_number = as.numeric(str_extract(Well, "\\d+"))
  ) %>%
  arrange(well_number, well_letter) %>%
  select(-well_letter, -well_number) %>% arrange(cell_line)

replicate_summary <- w1_cell_area %>%
  mutate(
    #mM_NaBu = factor(mM_NaBu, levels = c(0, 0.5, 1, 5)),
    replicate = Well
  ) %>%
  group_by(cell_line, mM_NaBu, Well) %>%
  summarise(
    mean_area = mean(cell_area, na.rm = TRUE),
    #sd_area = sd(cell_area, na.rm = TRUE),
    #n_cells = n(),
    #se = sd_area / sqrt(n_cells),
    #ci = qt(0.975, df = n_cells - 1) * se,
    .groups = "drop"
  )

summary_df <- replicate_summary %>%
  group_by(cell_line, mM_NaBu) %>%
  summarise(
    mean = mean(mean_area),
    sd = sd(mean_area),
    n = n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

replicate_summary %>% ggplot(aes(x = mM_NaBu, y = mean_area))+geom_point()+facet_wrap(~cell_line)


ggplot(summary_df,
       aes(x = mM_NaBu, y = mean, group = cell_line)) +
  
  geom_line(linewidth = 0.8) +
  
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.1,
    linewidth = 0.7
  ) +
  
  geom_point(
    shape = 22,        # filled square
    size = 4,
    fill = "black",
    color = "black"
  ) +
  
  geom_jitter(
    data = replicate_summary,
    aes(x = mM_NaBu, y = mean_area),
    width = 0.05,
    height = 0,
    shape = 21,        # empty circle
    size = 2.5,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~cell_line) +
  
  labs(
    x = "NaBu (mM)",
    y = "Mean cell area"
  ) +
  
  theme_classic(base_size = 14) + ggtitle("Mean cell area as a function of mM NaBu")


## w2
plate1 <- read_csv(file.path(exp_cell_path, "exp1_w2_phase-contrast_area/plate1/combined_cell_area.csv"))
plate2 <- read_csv(file.path(exp_cell_path, "exp1_w2_phase-contrast_area/plate2/combined_cell_area.csv"))
combined_w2 <- bind_rows("plate1" = plate1, "plate2" = plate2, .id = "plate")
combined_w2 %<>% mutate(cell_area = inch_to_um*area/ncell, mM_NaBu = factor(mM_NaBu))
combined_w2 %>% ggplot(aes(x = uM_Resveratrol, y = cell_area, colour = mM_NaBu))+geom_point()
combined_w2_summary <- combined_w2 %>%
  group_by(uM_Resveratrol, mM_NaBu) %>%
  summarise(
    mean = mean(cell_area),
    sd = sd(cell_area),
    n = n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

w2_summary_2 <- combined_w2 %>%
  group_by(plate, well, mM_NaBu, uM_Resveratrol) %>%
  summarise(
    mean = mean(cell_area),
    sd = sd(cell_area),
    n = n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

w2_summary_2 <- w2_summary_2 %>%
  mutate(
    well_letter = str_extract(well, "^[A-Z]"),
    well_number = as.numeric(str_extract(well, "\\d+"))
  ) %>%
  arrange(well_number, well_letter) %>%
  select(-well_letter, -well_number) %>% arrange(plate)

ggplot(combined_w2_summary,
       aes(x = uM_Resveratrol, y = mean, col = mM_NaBu)) +
  
  geom_line(linewidth = 0.8) +
  
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 1,
    linewidth = 0.8
  ) +
  
  geom_point(aes(fill=mM_NaBu),
    shape = 21,        # filled square
    size = 3,
    #fill = mM_NaBu,
    color = "black"
  ) +
  
  geom_jitter(
    data = combined_w2,
    aes(x = uM_Resveratrol, y = cell_area, color = mM_NaBu),
    width = 0.1,
    height = 0,
    shape = 21,        # empty circle
    size = 1.25,
    alpha = 0.30,
    #fill = "white",
    #color = "black",
    inherit.aes = FALSE
  ) +
  labs(
    #x = "NaBu (mM)",
    y = "Mean cell area"
  ) + 
  theme_classic(base_size = 14) +
  #theme(legend.position = "top") + 
  ggtitle("Mean cell area as a function of µM Resveratrol for DLD-1")
### coulter counter

coulter <- read_csv(file.path(exp_cell_path, "exp1_w1_dose_response_coulter/multisizer_summary.csv"))


dose_map <- c(
  "1" = 0,
  "2" = 0.5,
  "3" = 1,
  "4" = 5
)

coulter <- coulter %>%
  mutate(
    file_clean = str_remove(basename(file), "\\.[^.]+$"),
    file_clean = str_replace_all(file_clean, "\\s+", ""),
    
    well = str_extract(file_clean, "^[A-Za-z]\\d+|^\\d+[A-Za-z]"),
    
    replicate = case_when(
      str_detect(well, "^[A-Za-z]\\d+$") ~ str_extract(well, "^[A-Za-z]"),
      str_detect(well, "^\\d+[A-Za-z]$") ~ str_extract(well, "[A-Za-z]$"),
      TRUE ~ NA_character_
    ),
    
    well_number = case_when(
      str_detect(well, "^[A-Za-z]\\d+$") ~ str_extract(well, "\\d+"),
      str_detect(well, "^\\d+[A-Za-z]$") ~ str_extract(well, "^\\d+"),
      TRUE ~ NA_character_
    ),
    
    cell_line = str_split_i(file_clean, "_", 2),
    
    nabu_concentration = unname(
      as.numeric(dose_map[well_number])
    )
  ) %>%
  select(-file_clean, -well)

coulter %<>% mutate(estimated_area=pi*(mean_diameter/2)^2, 
                    cell_line = ifelse(str_detect(cell_line, "DLD"), "DLD1", cell_line)) %>% 
  mutate(cell_line = ifelse(str_detect(cell_line, "Caco-2"), "caco2", cell_line))
coulter %<>% rename("mM_NaBu" = nabu_concentration)
coulter %>% ggplot(aes(x = mM_NaBu, y = estimated_area, colour = cell_line))+geom_point()

coulter %>% ggplot(aes(x = mM_NaBu, y = cell_count, colour = cell_line))+geom_point()+ggtitle("coulter counter counts vs mM NaBu")

data %<>% mutate(well_cell = paste0(well_row, well_column, "_", cell_line))
coulter %<>% mutate(well_cell = paste0(replicate, well_number, "_", cell_line))

well_means <- data %>%
  group_by(
    well_row,
    well_column,
    cell_line
  ) %>%
  summarise(
    mean_area = mean(Area, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>% mutate(well_cell = paste0(well_row, well_column, "_", cell_line))

joined_data <- coulter %>% select(well_cell, estimated_area) %>% inner_join(well_means, by = "well_cell")

joined_data %<>% rename("coulter_area" = estimated_area, "imagej_area" = mean_area)

joined_data %>% ggplot(aes(x=coulter_area, y = imagej_area, colour = cell_line))+
  geom_point() + geom_smooth() + ggtitle("imagej and coulter area correlation")


### proliferation assay
prol <- read_csv(file.path(exp_cell_path, "proliferation_assay.csv")) %>% rename(RLU=Value, treatment=Condition)

summary_df <- prol %>%
  group_by(treatment, Day) %>%
  summarise(
    mean = mean(RLU),
    sd = sd(RLU),
    .groups = "drop"
  )

# Plot
ggplot(summary_df, aes(x = Day, y = mean, color = treatment, group = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2,
                linewidth = 0.7) +
  scale_x_continuous(breaks = 1:8) +
  labs(
    x = "Day",
    y = "Cell count",
    color = "treatment"
  ) +
  theme_classic(base_size = 14)
prol %>% ggplot(aes(x=Day, y=RLU, col=treatment))+geom_point()


caco_luciferase <- read_csv(file.path(exp_cell_path, "luciferase_assay_caco2.csv"))
caco_luciferase %<>% rename(treatment = factor("Condition"), RLU = "Value")
caco_luciferase %<>% mutate(replicate = gsub("[a-z]", "", Replicate, ignore.case=T))
caco_luciferase %>% ggplot(aes(x=treatment, y = RLU, col = replicate ))+geom_jitter(width=0.1)


# Calculate summary statistics
summary_caco_luciferase <- caco_luciferase %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(RLU),
    sd = sd(RLU),
    .groups = "drop"
  )

# Plot
ggplot(summary_caco_luciferase, aes(x = treatment, y = mean)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2,
                linewidth = 0.7) +
  geom_jitter(data = caco_luciferase,
              aes(x = treatment, y = RLU, fill = replicate),
              width = 0.08,
              size = 2.8,
              shape = 21,
              #fill = "replicate",
              color = "black",
              inherit.aes = FALSE) +
  labs(
    x = NULL,
    y = "Relative Light Unit"
  ) +
  theme_classic(base_size = 14) + 
  geom_text(
    aes(
      y = mean + sd,
      label = scales::label_scientific(digits = 3)(mean)
    ),
    vjust = -0.8,
    fontface = "bold",
    size = 4
  ) + ggtitle("Caco-2 luciferase assay")



### exp 4
w1_nuclei <- read_csv(file.path(exp_cell_path, "exp4_w1_nuclei.csv"))
w1_nuclei %<>% mutate(condition = factor(condition, levels = c("control", "NaBu", "NaBu_Resv")))
w1_nuclei %>% ggplot(aes(x=condition, y = nucleus_area))+geom_point()

w1_nuclei %<>% mutate(replicate = factor(if_else(image %% 2 == 1, 1, 2)))

summary_w1_nuclei <- w1_nuclei %>%
  group_by(condition) %>%
  summarise(
    mean = mean(nucleus_area),
    sd = sd(nucleus_area),
    .groups = "drop"
  )

ggplot(summary_w1_nuclei, aes(x = condition, y = mean)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.2,
    linewidth = 0.7
  ) +
  geom_jitter(
    data = w1_nuclei,
    aes(x = condition, y = nucleus_area, fill = replicate),
    width = 0.1,
    size = 2.5,
    shape = 21,
    #fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = round(mean, 2), y = mean + sd),
    vjust = -0.8,
    fontface = "bold",
    size = 4
  ) +
  labs(
    x = NULL,
    y = expression("Nucleus area ("*mu*m^2*")"),
  ) +
  theme_classic(base_size = 14) +
  #theme(legend.position = "none") +
  expand_limits(y = max(summary_w1_nuclei$mean + summary_w1_nuclei$sd) * 1.15) +ggtitle("Live imaging nucleus quantification")

# w2
w2_nuclei <- read_csv(file.path(exp_cell_path, "exp4_w2_nuclei.csv"))
w2_nuclei %<>% mutate(condition = factor(condition, levels = c("control", "NaBu", "Resv", "NaBu_Resv")))
w2_nuclei %>% ggplot(aes(x=condition, y = nucleus_area))+geom_point()

summary_w2_nuclei <- w2_nuclei %>%
  group_by(condition) %>%
  summarise(
    mean = mean(nucleus_area),
    sd = sd(nucleus_area),
    .groups = "drop"
  )

ggplot(summary_w2_nuclei, aes(x = condition, y = mean)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.2,
    linewidth = 0.7
  ) +
  geom_jitter(
    data = w2_nuclei,
    aes(x = condition, y = nucleus_area),
    width = 0.05,
    size = 2.5,
    shape = 21,
    #fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = round(mean, 2), y = mean + sd+ 20),
    vjust = -0.8,
    fontface = "bold",
    size = 4
  ) +
  labs(
    x = NULL,
    y = expression("Nucleus area ("*mu*m^2*")"),
  ) +
  theme_classic(base_size = 14) +
  #theme(legend.position = "none") +
  expand_limits(y = max(summary_w2_nuclei$mean + summary_w2_nuclei$sd) * 1.15) +ggtitle("Fixed cells nucleus quantification")


#### experiment 2
resazurin <- read_csv(file.path(exp_cell_path, "exp2_w2_b2_resazurin.csv"))
hoechst <- read_csv(file.path(exp_cell_path, "exp2_w2_b2_hoechst.csv"))
combined_fluo <- bind_rows("resazurin" = resazurin,"hoechst" = hoechst, .id = "fluorophore")
treatment_concentrations <- setNames(c("NaBu_1mM", "Resv_40µM", "NaBu_1mM_Resv_40µM"),c("NaBu", "Resv", "NaBu_Resv"))
#combined_fluo %<>% mutate(condition_old = condition, condition = ifelse(condition %in% names(treatment_concentrations), treatment_concentrations[condition], condition))
#combined_fluo %<>% mutate(condition = factor(condition, levels = c("media", "control", "NaBu_1mM", "Resv_40µM", "NaBu_1mM_Resv_40µM")))
combined_fluo %<>% mutate(condition = factor(condition, levels = c("media", "control", "NaBu", "Resv", "NaBu_Resv")))
combined_fluo %<>% filter(!(column %in% c(2,11)))

media_means <- combined_fluo %>%
  filter(condition == "media") %>%
  group_by(fluorophore, cell_line) %>%
  summarise(
    media_mean = mean(value),
    .groups = "drop"
  )
corrected_fluo <- combined_fluo %>%
  left_join(media_means,
            by = c("fluorophore", "cell_line")) %>%
  mutate(
    value_corrected = value - media_mean
  )

control_means <- corrected_fluo %>%
  filter(condition == "control") %>%
  group_by(fluorophore, cell_line) %>%
  summarise(
    control_mean = mean(value_corrected),
    .groups = "drop"
  )

fluo_normalized <- corrected_fluo %>%
  left_join(control_means,
            by = c("fluorophore", "cell_line")) %>%
  mutate(
    relative_intensity = value_corrected / control_mean
  )

# Summary statistics
summary_df <- fluo_normalized %>%
  group_by(cell_line, fluorophore, condition) %>%
  summarise(
    mean = mean(relative_intensity),
    sd = sd(relative_intensity),
    n = n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

ggplot(summary_df %>% filter(condition != "media"),
       aes(x = condition,
           y = mean,
           fill = fluorophore)) +
  
  geom_col(position = position_dodge(0.8),
           width = 0.7,
           colour = "black") +
  
  geom_errorbar(
    aes(ymin = mean - sem,
        ymax = mean + sem),
    width = 0.2,
    position = position_dodge(0.8)
  ) +
  
  geom_jitter(
    data = fluo_normalized %>% filter(condition != "media"),
    aes(x = condition,
        y = relative_intensity,
        fill = fluorophore),
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.8
    ),
    shape = 21,
    size = 2,
    colour = "black",
    inherit.aes = FALSE
  ) + geom_text(
    aes(
      label = round(mean, 2),
      y = mean + 0.85 # + sem
    ),
    position = position_dodge(width = 0.8),
    vjust = -0.8,
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(~cell_line) +
  
  labs(
    x = NULL,
    y = "Relative Fluorescence"
  ) +
  
  theme_classic(base_size = 14) + theme(
    #axis.text.x = element_text(angle = -45, hjust = 1),
    legend.position = "top"
  ) + ggtitle("Effect on relative cell viability upon treatment with NaBu and/or Resv")

# ## old plot
# summary_df_old <- combined_fluo %>%
#   group_by(cell_line, fluorophore, condition) %>%
#   summarise(
#     mean = mean(value),
#     sd = sd(value),
#     n = n(),
#     sem = sd / sqrt(n),
#     .groups = "drop"
#   )
# 
# ggplot(summary_df_old,
#        aes(x = condition,
#            y = mean,
#            fill = fluorophore)) +
#   
#   geom_col(position = position_dodge(0.8),
#            width = 0.7,
#            colour = "black") +
#   
#   geom_errorbar(
#     aes(ymin = mean - sem,
#         ymax = mean + sem),
#     width = 0.2,
#     position = position_dodge(0.8)
#   ) +
#   
#   geom_jitter(
#     data = combined_fluo,
#     aes(x = condition,
#         y = value,
#         fill = fluorophore),
#     position = position_jitterdodge(
#       jitter.width = 0.1,
#       dodge.width = 0.8
#     ),
#     shape = 21,
#     size = 2,
#     colour = "black",
#     inherit.aes = FALSE
#   ) +
#   
#   facet_wrap(~cell_line, scales = "free_y") +
#   
#   labs(
#     x = NULL,
#     y = "Relative Fluorescence"
#   ) +
#   
#   theme_classic(base_size = 14)


