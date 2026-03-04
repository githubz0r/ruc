library(tidyverse)
library(magrittr)

resin_loading <- 0.24 #mmol/g

data <- read_excel('/Users/b246357/Documents/uni/peptides/ptb2026uvfmoc.xlsx')
data_clean_header <- data
data_clean_header %<>% rename('wavelength' = ...1); data_clean_header <- data_clean_header[2:nrow(data_clean_header), ]
data_clean_long <- data_clean_header %>% pivot_longer(cols=-wavelength, values_to = 'abs', names_to = 'group')

wl_301 <- data_clean_long %>% filter(wavelength == '301')

get_conc <- function(absorbance, mg_beads){
  conc_1 <- 10*absorbance/(7800) # was diluted in dmf 10x, mol/l
  total_fmoc <- conc_1 *3*10^(-3) # we took 1 ml from 3 ml, mol
  mmol_g_resin <- 10^3 * total_fmoc/(mg_beads*10^(-3))
  return(mmol_g_resin)
}
#get_conc(0.438-8.9999999999999993E-3)
g1_7.9 <- get_conc(0.438 - 8.9999999999999993E-3, 7.9)
g1_13.5 <- get_conc(0.81599999999999995 - 8.9999999999999993E-3, 13.5)

g3_14 <- get_conc(0.69299999999999995 - 3.0000000000000001E-3, 14)
g3_15.8 <- get_conc(0.80300000000000005 - 3.0000000000000001E-3, 15.8)

g4_14 <- get_conc(0.83799999999999997 - 2E-3, 14)
g4_9.8 <- get_conc(0.55900000000000005 - 2E-3, 9.8)

measurements <- data.frame(g1=c(g1_7.9, g1_13.5), g3=c(g3_14, g3_15.8), g4=c(g4_14, g4_9.8)) %>% pivot_longer(cols=everything(), values_to = 'conc', names_to = 'group')
measurements$conc %>% summary()
                           