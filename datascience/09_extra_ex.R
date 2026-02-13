library(tidyverse)
library(ggplot2)
library(magrittr)
library(emmeans)
library(multcomp)
library(multcompView)

wisteria <- c("grey65", "burlywood3", "khaki2", "plum1", "lightcyan2", "cornflowerblue", "slateblue3")
chicken_experiment <- read_csv("https://raw.githubusercontent.com/johnshorter/TeachingData/main/betacarotenoid.csv")

ggplot(chicken_experiment, aes(x = Pasture, y = Beta)) +
  # 1. Boxplot in the background (hide outliers to avoid doubling)
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  #facet_wrap(~Field, scales = "free_y") +
  
  # 2. Individual points on top (jittered)
  geom_jitter(aes(fill = Field), shape=21,  width = 0.1, alpha = 0.5, color = "black") +
  #geom_point(aes(fill = treatment), shape = 21, color = "black", alpha = 0.8,
  #position = position_jitter(width = 0.1, seed = 666)) + 
  
  # 3. Apply the palette
  #scale_fill_manual(values = wisteria) +
  
  # 4. Clean up the axis and theme
  theme_minimal() +
  theme(
    axis.title.x = element_blank(), # Remove "treatment" label
    legend.position = "bottom"        # Hide legend (x-axis labels are enough)
  ) +
  labs(
    title = "Distribution of Beta by Pasture",
    y = "Beta"
    #x = "Pasture"
  )

# if we had more values and some of them were numerical factors we could do lines instead
# of coloring by category or doing a facet_wrap
