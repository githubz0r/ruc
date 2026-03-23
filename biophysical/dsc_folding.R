library(tidyverse)
library(magrittr)
library(readxl)
library(broom)
library(fitdistrplus)
barnase <- read_excel("/Users/b246357/Documents/uni/biophysical/barnase_ubiquitin.xlsx", sheet = "Barnase")
ubiquitin <- read_excel("/Users/b246357/Documents/uni/biophysical/barnase_ubiquitin.xlsx", sheet = "Ubiquitin")
combined <- bind_rows(barnase=barnase, ubiquitin=ubiquitin, .id='protein')
colnames(combined) <- c("protein", "kelvin", "cp_kj_Kmol")

combined %>% ggplot(aes(x=kelvin, y=cp_kj_Kmol, colour = protein))+geom_point()+scale_x_continuous(n.breaks = 10)

barnase_cp_n <- combined %>%  filter(protein == 'barnase' & kelvin < 290)
barnase_cp_u <- combined %>%  filter(protein == 'barnase' & kelvin > 319)
barnase_cp_df <- bind_rows(barnase_cp_n, barnase_cp_u) %>% mutate(norm_k = kelvin - min(kelvin))
ub_cp_n <- combined %>%  filter(protein == 'ubiquitin' & kelvin < 326)
ub_cp_u <- combined %>%  filter(protein == 'ubiquitin' & kelvin > 370)
ub_cp_df <- bind_rows(ub_cp_n, ub_cp_u) %>% mutate(norm_k = kelvin - min(kelvin))

baselines <- list(barnase_cp_n = barnase_cp_n, 
                  barnase_cp_u = barnase_cp_u, 
                  ub_cp_n = ub_cp_n,
                  ub_cp_u = ub_cp_u)
cp_fits <- baselines %>%
  map(~ lm(cp_kj_Kmol ~ kelvin, data = .x)) %>%  # Run model for each DF
  map_dfr(tidy, .id = "dataset_name")
cps <- c(barnase_n = 17.5, barnase_u = 25, ub_n = 10, barnase_u = 20)


combined %<>% mutate(excess_cp = ifelse(protein == 'barnase', cp_kj_Kmol - cps['barnase_n'], cp_kj_Kmol - cps['ub_n']))
#combined_long <- combined %>% pivot_longer(cols = c(cp_kj_Kmol, excess_cp), values_to = 'cp', names_to = 'cp_type') 

#combined_long %>% ggplot(aes(x=kelvin, y=cp, colour = protein))+geom_point()+
#  scale_x_continuous(n.breaks = 10) + 
#  facet_grid(~cp_type)


# a: lower asymptote, d: upper asymptote, c: inflection point, b: hill slope

barnase_sigmoid <- nls(
  cp_kj_Kmol ~ a + (d - a) / (1 + exp(b * (c - norm_k))),
  data = barnase_cp_df,
  start = list(a = 20, d = 25, b = 0.1, c = 40) # Initial guesses
)


ubiquitin_sigmoid <- nls(
  cp_kj_Kmol ~ a + (d - a) / (1 + exp(b * (c - norm_k))),
  data = ub_cp_df,
  start = list(a = 10, d = 20, b = 0.05, c = 50) # Initial guesses
  ,nls.control(maxiter = 1000))

pred_barnase_sigmoid <- barnase_sigmoid %>% predict(newdata = barnase %>% mutate(norm_k = `T (K)` - min(barnase_cp_df$kelvin)))
pred_ub_sigmoid <- ubiquitin_sigmoid %>% predict(newdata = ubiquitin %>% mutate(norm_k = `T (K)` - min(ub_cp_df$kelvin)))

#pred_df <- bind_rows(data.frame(protein='barnase', cp_sigmoid = pred_barnase_sigmoid),
#                     data.frame(protein='ubiquitin', cp_sigmoid = pred_ub_sigmoid))

combined %<>% mutate("sigmoid_cp" = c(pred_barnase_sigmoid, pred_ub_sigmoid))
combined_long2 <- combined %>% pivot_longer(cols = c(cp_kj_Kmol, sigmoid_cp), values_to = 'cp', names_to = 'cp_type') 
combined_long2 %>% ggplot(aes(x=kelvin, y=cp, colour = protein, shape = cp_type))+geom_point()+
  scale_x_continuous(n.breaks = 10) #+ 
  #facet_grid(~cp_type)

combined %<>% mutate(corrected_cp = cp_kj_Kmol - sigmoid_cp)
combined %>% ggplot(aes(x = kelvin, y = corrected_cp, colour = protein))+geom_point()

### testing numerical integration
a <- 300
b <- 330

number_of_rectangles <- 44
xgrid <- seq(a, b, length.out=number_of_rectangles+1)
dx <- (b-a)/number_of_rectangles
#xs <- xgrid[1:number_of_rectangles]
ys <- combined %>% filter(protein == 'barnase' & (kelvin >= 300 & kelvin <= 330)) %>% pluck('corrected_cp')
sum_areas <- dx*sum(ys)
###

gaussian_fun <- function(x, k, mu, sigma) {
  k * exp(-((x - mu)^2) / (2 * sigma^2))
}

barnase_fit_norm <- nls(corrected_cp ~ gaussian_fun(kelvin, k, mu, sigma), 
                  data = combined %>% filter(protein == 'barnase'),
                  start = list(k = 45, 
                             mu = 300, 
                             sigma = 10))
ubiquitin_fit_norm <- nls(corrected_cp ~ gaussian_fun(kelvin, k, mu, sigma), 
                        data = combined %>% filter(protein == 'ubiquitin'),
                        start = list(k = 15, 
                                     mu = 340, 
                                     sigma = 10))

barnase_params <- coef(barnase_fit_norm)
ubiquitin_params <- coef(ubiquitin_fit_norm)

barnase_fitted <- function(x) {
  gaussian_fun(x, k = barnase_params["k"], mu = barnase_params["mu"], sigma = barnase_params["sigma"])
}
ubiquitin_fitted <- function(x) {
  gaussian_fun(x, k = ubiquitin_params["k"], mu = ubiquitin_params["mu"], sigma = ubiquitin_params["sigma"])
}

barnase_gaussian_pred <- barnase_fitted(combined %>% filter(protein == 'barnase') %>% pluck('kelvin'))
ubiquitin_gaussian_pred <- ubiquitin_fitted(combined %>% filter(protein == 'ubiquitin') %>% pluck('kelvin'))

combined %<>% mutate(gaussian_fits = c(barnase_gaussian_pred, ubiquitin_gaussian_pred))

combined_long3 <- combined %>% pivot_longer(cols = c(corrected_cp, gaussian_fits), values_to = 'corrected_cp', names_to = 'cp_type') %>% 
  mutate('cp_or_pred_cp' = paste0(protein, '_', cp_type))
combined_long3 %>% ggplot(aes(x=kelvin, y=corrected_cp, colour = cp_or_pred_cp))+geom_point()+geom_line()
  scale_x_continuous(n.breaks = 10)

d_h_cal_barnase <- integrate(barnase_fitted, lower = 200, upper = 400)
d_h_cal_ubiquitin <- integrate(ubiquitin_fitted, lower = 200, upper = 400)

#tm_barnase <- combined %>% filter(protein == 'barnase') %>% filter(cp_kj_Kmol == max(cp_kj_Kmol)) %>% pluck('kelvin')
#tm_ubiquitin <- combined %>% filter(protein == 'ubiquitin') %>% filter(cp_kj_Kmol == max(cp_kj_Kmol)) %>% pluck('kelvin') 
# I don't need these i have gaussians

# guessing delta cp based on values from sigmoid fits
d_cp_barnase <- coef(barnase_sigmoid)['d'] - coef(barnase_sigmoid)['a']
d_cp_ubiquitin <- coef(ubiquitin_sigmoid)['d'] - coef(ubiquitin_sigmoid)['a']

mt_barnase <- barnase_params['mu']
mt_ubiquitin <- ubiquitin_params['mu']



d_s_tm_barnase <- d_h_cal_barnase$value/barnase_params['mu'] # -1.1 kj/mol*K
d_s_tm_ubiquitin <- d_h_cal_ubiquitin$value/ubiquitin_params['mu'] # -0.62 kj/mol*K

# function corresponding to 27.28 in klostermeier
enthalpy <- function(kelvin, m_enthalpy, mt, d_cp){
  d_h <- m_enthalpy + d_cp*(kelvin - mt)
  return(d_h)
}

# 27.29 in klostermeier
delta_g <- function(kelvin, m_enthalpy, mt, d_cp){
  d_g <- m_enthalpy * (1 - kelvin/mt) + d_cp * (kelvin - mt - kelvin * log(kelvin/mt))
  return(d_g)
  #min_t_s <- d_g - enthalpy(kelvin, m_enthalpy, mt, d_cp)
}

entropy <- function(kelvin, m_entropy, mt, d_cp){
  d.s <- m_entropy + d_cp*(log(kelvin/mt))
  return(d.s)
}
thermo_params <- function(kelvin, m_enthalpy, m_entropy, mt, d_cp){
  d_g <- delta_g(kelvin, m_enthalpy, mt, d_cp)
  d_h <- enthalpy(kelvin, m_enthalpy, mt, d_cp)
  t_d_s <- -(d_g - d_h)
  d_s2 <- entropy(kelvin, m_entropy, mt, d_cp)
  d_g2 <- d_h - kelvin*d_s2
  return(list(d_g = d_g, d_h = d_h, t_d_s = t_d_s, t_d_s2 = kelvin*d_s2, d_g2 = d_g2))
}

kelvin_range <- 290:390
barnase_thermos <- kelvin_range %>% 
  lapply(thermo_params, d_h_cal_barnase$value, d_s_tm_barnase, mt_barnase, d_cp_barnase) %>% 
  bind_rows() %>% mutate(kelvin = kelvin_range)

ubiquitin_thermos <- kelvin_range %>% 
  lapply(thermo_params, d_h_cal_ubiquitin$value, d_s_tm_ubiquitin, mt_ubiquitin, d_cp_ubiquitin) %>% 
  bind_rows() %>% mutate(kelvin = kelvin_range)

thermos_combined <- bind_rows(barnase = barnase_thermos, ubiquitin = ubiquitin_thermos, .id = 'protein')
thermos_combined_long <- thermos_combined %>% pivot_longer(cols=c(d_g, d_h, t_d_s), values_to = 'thermval', names_to = 'thermparam')


# keep in mind this is delta G for the unfolded state
thermos_combined_long %>% ggplot(aes(y=thermval, x = kelvin, col = thermparam, shape = protein)) + 
  geom_point()

thermos_combined %>% ggplot(aes(x=kelvin, y=d_g, colour=protein))+geom_point()

### van't hoff analysis
# loop over temperature values and integrate
fu_t <- function(kelvin, gaussian_function, full.integral, lower = 200){
  integral <- integrate(gaussian_function, lower = lower, upper = kelvin) %>% pluck('value')
  f.u <- integral/full.integral
  return(f.u)
  #Ku <- f.u/(1-f.u)
}

barnase_fu_t <- combined %>% filter(protein == 'barnase') %>% pluck('kelvin') %>% 
  sapply(fu_t, barnase_fitted, d_h_cal_barnase$value)

ubiquitin_fu_t <- combined %>% filter(protein == 'ubiquitin') %>% pluck('kelvin') %>% 
  sapply(fu_t, ubiquitin_fitted, d_h_cal_ubiquitin$value)

combined %<>% mutate(f_u = c(barnase_fu_t, ubiquitin_fu_t))

combined %>% ggplot(aes(x=kelvin, y=f_u, colour=protein)) + geom_point()
combined <- combined %>%
  group_by(protein) %>%
  mutate(relative_kelvin = kelvin - min(kelvin, na.rm = TRUE)) %>%
  ungroup()

combined %<>% mutate(ln_ku = log(f_u/(1-f_u)))
combined_no_nan <- combined %>% filter(!is.na(ln_ku)) %>% mutate(inv_kelvin = 1/kelvin) # a few values became nan by the transformation
combined_no_nan %<>% filter(!(protein == "barnase" & kelvin > 330)) # some values mess with the fit because of floating point imprecision

#combined_no_nan %<>% filter(inv_kelvin > 0.00275 & inv_kelvin < 0.003)

vant_hoff_barnase <- lm(ln_ku ~ inv_kelvin, data = combined_no_nan %>% filter(protein == 'barnase'))
vant_hoff_barnase_entropy <- coef(vant_hoff_barnase)['(Intercept)'] * 8.314 * 10^-3
vant_hoff_barnase_enthalpy <- -coef(vant_hoff_barnase)['inv_kelvin'] * 8.314 * 10^-3
pred_vant_hoff_barnase <- predict(vant_hoff_barnase)
vant_hoff_ubiquitin <- lm(ln_ku ~ inv_kelvin, data = combined_no_nan %>% filter(protein == 'ubiquitin'))
vant_hoff_ubiquitin_entropy <- coef(vant_hoff_ubiquitin)['(Intercept)'] * 8.314 * 10^-3
vant_hoff_ubiquitin_enthalpy <- -coef(vant_hoff_ubiquitin)['inv_kelvin'] *8.314 * 10^-3
pred_vant_hoff_ubiquitin <- predict(vant_hoff_ubiquitin)
combined_no_nan %<>% mutate(pred_vant_hoff = c(pred_vant_hoff_barnase, pred_vant_hoff_ubiquitin))

combined_no_nan_long <- combined_no_nan %>% pivot_longer(cols = c(ln_ku, pred_vant_hoff), values_to = "lnK", names_to = "emp_or_calc")

combined_no_nan_long %>% ggplot(aes(x=inv_kelvin, y = lnK, fill = protein, colour = emp_or_calc, shape = emp_or_calc))+
  geom_point() + 
  scale_shape_manual(values = c(21, 22)) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,    # Use a simple circle for the protein legend
        color = NA,    # REMOVE the black border so the fill is visible
        stroke = 0     # Set border thickness to zero
      )
    ),
    color = guide_legend(
      override.aes = list(
        fill = "white" # Set a neutral background for the 'type' shapes
      )
    ),
    shape = guide_legend() # Keep this linked to color to merge them
  ) +
  theme_minimal() + 
  labs(title = "vant hoff",
       fill = "Protein Name",
       color = "Measurement Type",
       shape = "Measurement Type")
# try nonlinear instead

# vant_hoff_exp_function <- function(kelvin, enthalpy, entropy){
#   K <- exp(-enthalpy/(kelvin*8.314) + entropy/8.314)
#   return(K)
# }
# 
# combined %<>% mutate(ku = f_u/(1-f_u))
# barnase_vant_hoff_exp <- nls(ku ~ vant_hoff_exp_function(kelvin, enthalpy, entropy), 
#                         data = combined %>% filter(protein == 'barnase'),
#                         start = list(enthalpy = 262,
#                                      entropy = 1.1))

## barnase
# barnase variables dsc
barnase_thermos_298 <- thermos_combined %>% filter(protein == 'barnase' & kelvin == 298) %>% mutate(d_s = t_d_s/298)
print(barnase_thermos_298)
# barnase variables vant hoff
print(vant_hoff_barnase_enthalpy %>% unname()) 
print(vant_hoff_barnase_entropy %>% unname())
print((vant_hoff_barnase_enthalpy - 298*vant_hoff_barnase_entropy) %>% unname())


# ubiquitin variables vant hoff
ubiquitin_thermos_298 <- thermos_combined %>% filter(protein == 'ubiquitin' & kelvin == 298) %>% mutate(d_s = t_d_s/298)
print(ubiquitin_thermos_298)
print(vant_hoff_ubiquitin_enthalpy %>% unname()) 
print(vant_hoff_ubiquitin_entropy %>% unname())
print((vant_hoff_ubiquitin_enthalpy - 298*vant_hoff_ubiquitin_entropy) %>% unname())
  
