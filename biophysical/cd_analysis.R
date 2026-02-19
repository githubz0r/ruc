library(tidyverse)
library(magrittr)
nap <- read_delim('/Users/b246357/Documents/uni/biophysical/grp1_nap_protein_min_blank.csv', delim = ',', skip = 19)
### remove the headers and tail values
nap <- nap[1:131, ]
nap_clean <- nap %>%
  separate(1, into = c("wavelength", "theta", "ht"), sep = ",") %>%
  mutate(across(everything(), as.numeric))

acetate <- read_delim('/Users/b246357/Documents/uni/biophysical/grp1_acetate_protein_mink_blank.csv', delim = ',', skip = 19)
acetate <- acetate[1:131, ]
acetate_clean <- acetate %>%
  separate(1, into = c("wavelength", "theta", "ht"), sep = ",") %>%
  mutate(across(everything(), as.numeric))

mw <- 14313 # g/mol
conc_nap <- 0.091 # mg/ml
conc_acetate <- 0.085
l_ <- 0.2 # cm
naa <- 129
conc_micromolar_nap <- 10^6 * conc_nap/mw 
conc_micromolar_acetate <- 10^6 * conc_acetate/mw

mean_resid_theta <- function(theta, concentration){
  return(theta * mw / (10 * naa * concentration * l_))
}
nap_clean %<>% mutate(mrt = mean_resid_theta(theta, conc_nap))
acetate_clean %<>% mutate(mrt = mean_resid_theta(theta, conc_acetate))

nap_acetate <- bind_rows(nap = nap_clean, acetate = acetate_clean, .id = 'buffer')

nap_acetate %>% ggplot(aes(x=wavelength, y = theta, color = buffer))+geom_point()
nap_acetate %>% ggplot(aes(x=wavelength, y = ht, color = buffer))+geom_point()

helix_fraction_nap <- (nap_clean %>% filter(wavelength == 222) %>% pluck('mrt'))/(-33000)
helix_fraction_acetate <- (acetate_clean %>% filter(wavelength == 222) %>% pluck('mrt'))/(-33000)
print(c(helix_fraction_acetate, helix_fraction_nap))