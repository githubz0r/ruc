library(tidyverse)
library(magrittr)
library(minpack.lm)

heats <- c(-41.84, -36.82, -23.01, -18.83, -11.72, -7.95, -7.12, -6.69, -3.35, -2.93, -1.67, -1.67, -1.26, -0.84, -0.42)
molar_ratios <- c(0.27394258, 0.54788516, 0.82182774, 1.09577033, 1.36971291, 1.64365549, 1.91759807, 2.19154065, 2.46548323, 2.73942582, 3.0133684, 3.28731098, 3.56125356, 3.83519614, 4.10913872)
volumes <- c(0.2053, 0.2078, 0.2103, 0.2128, 0.2153, 
                    0.2178, 0.2203, 0.2228, 0.2253, 0.2278, 
                    0.2303, 0.2328, 0.2353, 0.2378, 0.2403) * 10^-3

c_ligand <- c(
  0.012, 0.024, 0.036, 0.047, 0.058, 
  0.069, 0.079, 0.090, 0.100, 0.110, 
  0.119, 0.129, 0.138, 0.147, 0.156
) * 10^-3

c_protein <- c(
  0.04445, 0.04392, 0.04340, 0.04289, 0.04239, 
  0.04190, 0.04143, 0.04096, 0.04051, 0.04006, 
  0.03963, 0.03920, 0.03878, 0.03838, 0.03798
) * 10^-3

qi <- c(
  -1.06799e-07, -1.92729e-07, -2.56379e-07, -3.00949e-07, -3.31709e-07,
  -3.53284e-07, -3.68867e-07, -3.80486e-07, -3.89408e-07, -3.96438e-07,
  -4.02102e-07, -4.06752e-07, -4.10633e-07, -4.13917e-07, -4.16730e-07
)
itc_df <- data.frame(heat=heats, molar_ratio=molar_ratios, volumes=volumes, c_ligand=c_ligand, c_protein=c_protein, qi = qi)

protein_conc <- 0.045*(10^-3)

itc_function <- function(molar_ratio, volume, enthalpy, r.w, Nl){
  heat <- 0.5*enthalpy*volume*(
    1-(Nl*molar_ratio -1 + r.w)/sqrt(
      (Nl*molar_ratio)^2 - 2*Nl*molar_ratio*(1-r.w) + (1+r.w^2)
      )
    )
  return(heat)
}

itc_function2 <- function(c_prot, c_l, volume, enthalpy, kd, Nl){
  heat <- 0.5*enthalpy*volume*(1-(
    Nl*c_l - c_prot + kd)/sqrt(
      (Nl*c_l + c_prot + kd)^2 - 4*c_prot*Nl*c_l
      )
    )
  return(heat)
}

itc_function3 <- function(c_prot, c_l, volume, enthalpy, kd, Nl){
  qi <- enthalpy*volume*0.5*(Nl*c_l + c_prot + kd - sqrt(
      (Nl*c_l + c_prot + kd)^2 - 4*c_prot*Nl*c_l
    )
  )
  return(qi)
}


# itc_function_simple <- function(molar_ratio, enthalpy, r.w, Nl){
#   heat <- 0.5*enthalpy*0.0002053*(1-(Nl*molar_ratio -1 + r.w)/sqrt((Nl*molar_ratio)^2 - 2*Nl*molar_ratio*(1-r.w) + (1+r.w)^2))
#   return(heat)
# }
# 
# m_simple <- nlsLM(heat ~ itc_function_simple(molar_ratio, enthalpy, r.w, Nl), 
#            data = itc_df,
#            start = list(enthalpy = -50, r.w = 0.2, Nl = 1),
#            #lower = c(enthalpy = -100, kd = 10^(-9), Nl = 0.5), # Minimum possible values
#            #upper = c(enthalpy = -10, kd = 10^(-3), Nl = 2.0),
#            control = nls.lm.control(maxiter = 10000, 
#                                     maxfev = 20000, 
#                                     ftol = 1e-9))
# this dogshit function doesn't work, maybe the protein concentration needs to vaery?



m <- nlsLM(heat ~ itc_function(molar_ratio, volumes, enthalpy, r.w, Nl), 
         data = itc_df,
         start = list(enthalpy = -50, r.w = 0.02, Nl = 1),
         lower = c(enthalpy = -100, r.w = 0, Nl = 0.5), # Minimum possible values
         #upper = c(enthalpy = -10, kd = 10^(-3), Nl = 2.0),
         control = nls.lm.control(maxiter = 10000, 
                                  maxfev = 20000, 
                                  ftol = 1e-9)
         )

m_pred <- predict(m)

itc_df %<>% mutate(heat_pred = m_pred)
itc_df_long <- itc_df %>% pivot_longer(cols = c(heat, heat_pred), names_to = "exp_or_pred", values_to = 'heats')
itc_df_long %>% ggplot(aes(x=molar_ratio, y = heats, color = exp_or_pred))+geom_point()

r_w_exp <- 0.0000159/protein_conc
m_pred_2 <- itc_function(itc_df$molar_ratio, itc_df$volumes, -50.2, 0.0000159, 1.25)
#m_pred_2
itc_function2(itc_df$c_protein, itc_df$c_ligand, itc_df$volumes, -50.2, 0.0000159, 1.25)
itc_function3(itc_df$c_protein, itc_df$c_ligand, itc_df$volumes, -50.2, 0.0000159, 1.25)

m3 <- nlsLM(qi ~ itc_function3(c_protein, c_ligand, volumes, enthalpy, kd, Nl), 
           data = itc_df,
           start = list(enthalpy = -50, kd = 10^(-6), Nl = 1),
           lower = c(enthalpy = -100, kd = 0, Nl = 0.5), # Minimum possible values
           #upper = c(enthalpy = -10, kd = 10^(-3), Nl = 2.0),
           control = nls.lm.control(maxiter = 10000, 
                                    maxfev = 20000, 
                                    ftol = 1e-8)
)

# fraction bound: [PL]/[Ptotal], use either binding isotherm or the more complicated function to get it
