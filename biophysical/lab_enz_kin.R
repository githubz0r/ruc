library(tidyverse)
library(magrittr)
library(readxl)
concentration_letter_key <- setNames(c(1000, 500, 250, 125,62,31,15), LETTERS[1:7]) # conc in mM

nadh_e <- 6220 # M^-1 cm^-1
light_path = 0.78 # cm 
enzyme_conc = 0.1*10^-3 # mM

data <- read_excel("/Users/b246357/Documents/uni/biophysical/ex8data.xlsx", skip=2) #%>% rename(time_s = time)
data_g1 <- data %>% select(Time, matches("[A-G]{1}[2-3]{1}$"))
data_g2 <- data %>% select(Time, matches("[A-G]{1}[7-9]{1}$"))
data_g3 <- data %>% select(Time, matches("[A-G]{1}(9|10|11)\\b{1}$"))
data_g1 <- data_g1[1:7, ]
data_g2 <- data_g2[1:7, ]
data_g3 <- data_g3[1:7, ]
data_g1 %<>% mutate(Time = as.numeric(Time))
data_g2 %<>% mutate(Time = as.numeric(Time))
data_g3 %<>% mutate(Time = as.numeric(Time))
data_g1_long <- data_g1 %>% pivot_longer(cols=-c(Time), values_to = "absorbance", names_to = "well")
data_g2_long <- data_g2 %>% pivot_longer(cols=-c(Time), values_to = "absorbance", names_to = "well")
data_g3_long <- data_g3 %>% pivot_longer(cols=-c(Time), values_to = "absorbance", names_to = "well")

#data_small %>% ggplot(aes(y=A2, x=Time))+geom_point()
data_g1_long %>% ggplot(aes(y=absorbance, x=Time))+geom_point()+facet_wrap(.~well, nrow=2) # 7 looks good
data_g2_long %>% ggplot(aes(y=absorbance, x=Time))+geom_point()+facet_wrap(.~well, nrow=2) # 7 looks good
data_g3_long %>% ggplot(aes(y=absorbance, x=Time))+geom_point()+facet_wrap(.~well, nrow=2) # 7 looks good

regressions_g1 <- data_g1 %>% select(-Time) %>% lapply(function(absorbance){
  conc <- 10^3 * absorbance/(nadh_e * light_path) # convert absorbance to mM
  model <- lm(conc ~ data_g1$Time)
  coefs <- broom::tidy(model)
  return(coefs)
})


slopes_g1 <- regressions_g1 %>% lapply(function(x){
  slope <- x %>% filter((term!='(Intercept)')) %>% pluck('estimate')
  return(slope)
})

v_c_df <- slopes_g1 %>% names %>% lapply(function(x){
  conc <- concentration_letter_key[str_sub(x, 1, 1)]
  value <- slopes_g1[[x]]
  row.id <- str_sub(x, 2, 2)
  return(list(conc = conc, velocity = value, row.id = row.id))
}) %>% bind_rows()


#v_c_df %<>% mutate(corrected_velocity=velocity*10^3/(nadh_e*light_path))

v_c_df %>% ggplot(aes(x=conc, y=velocity, colour=row.id))+geom_point()
fit_nl <- nls(velocity ~ (Vmax * conc) / (Km + conc), 
           data = v_c_df, 
           start = list(Vmax = 10, Km = 50))

# View the results
summary(fit_nl)
coefs_nl <- coef(fit_nl)
predicted_nl <- bind_cols(velocity = predict(fit_nl), conc = v_c_df$conc, row.id = 'predicted')
v_c_df %>% bind_rows(predicted_nl) %>% ggplot(aes(x=conc, y=velocity, colour=row.id))+geom_point()

# check lineweaver burke
v_c_df %<>% mutate(inv_vel = 1/velocity, inv_conc = 1/conc)

lineweaver_burk <- lm(inv_vel ~ inv_conc, data = v_c_df)

summary(lineweaver_burk)
intercept <- coef(lineweaver_burk)[1]
slope <- coef(lineweaver_burk)[2]
v_max <- 1/intercept
km <- slope*v_max

predicted_linear <- bind_cols(inv_vel = predict(lineweaver_burk), inv_conc = v_c_df$inv_conc, row.id = 'lwburke')
v_c_df %>% bind_rows(predicted_linear) %>% ggplot(aes(x=inv_conc, y=inv_vel, colour=row.id))+geom_point()




