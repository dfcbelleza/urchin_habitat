# Data analysis for the manuscript:
# "The behavior of sympatric urchin species across and ecosystem state gradient"

# sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 11 (bullseye)


# Load packages for data reading, analysis, and plotting
library(tidyverse)  # For data wrangling
library(lubridate)  # For easy manipulation of dates
library(brms)       # For Bayesian models
library(tidybayes)  # For extracting tidy data from Bayesian model fits
library(bayesplot)  # For added plotting functionality of Bayesian models
library(showtext)   # For specifying text/graphics

########################### LOAD AND PREPARE DATA ##############################
df_urchin   = read_csv("./data/urchin_data.csv")         # Urchin size-weight model fit data
df_quad     = read_csv("./data/transect_quad_data.csv")  # Benthic quadrat data. Excludes "Others"
df_depth    = read_csv("./data/depth_data.csv")          # Depth logger data
df_light    = read_csv("./data/temp_ppfd.csv")           # Mean monthly temperature and light (ppfd) data
df_tag      = read_csv("./data/tag_urchin.csv")          # Mark-recapture experiment data
df_rugosity = read_csv("./data/raw_rugosity.csv")        # Benthic rugosity data
df_wave     = read_csv("./data/wave_data.csv")           # Wave logger data
urchin_orig = read_csv("./data/urchin_orig.csv")         # Monthly urchin monitoring data

# Some descriptive statistics -------------------------------------------------

## Observed average daily wave heights -----
wave = df_wave |> 
  mutate(month = month(datetime)) |> 
  mutate(season = ifelse(month <= 2 | month == 12, "winter",
                  ifelse(month >= 3 & month < 6, "spring", 
                  ifelse(month >= 6 & month < 9, "summer",
                  ifelse(month >= 9 & month < 12, "autumn", "none"))))) |> 
  mutate(season = factor(season, levels = c("winter", "spring", "summer", "autumn"))) |> 
  group_by(location, season) |> 
  summarise_if(is.numeric, mean, na.rm = TRUE) 

wave |> 
  filter(str_detect(location, "Iso")) |> 
  pull(mean_height) |> range()

wave |> 
  filter(str_detect(location, "Veg")) |> 
  pull(mean_height) |> range()


## Observed seasonal light and temperature
pendant = df_light |> 
  mutate(month = month(datetime)) |> 
  mutate(season = ifelse(month <= 2 | month == 12, "winter",
                  ifelse(month >= 3 & month < 6, "spring", 
                  ifelse(month >= 6 & month < 9, "summer",
                  ifelse(month >= 9 & month < 12, "autumn", "none"))))) 

# Averaged across locations
pendant |> 
  filter(str_detect(season, "wi")) |> 
  summarise_if(is.numeric, mean, na.rm = TRUE) |> 
  dplyr::select(-month)

pendant |> 
  filter(str_detect(season, "su")) |> 
  summarise_if(is.numeric, mean, na.rm = TRUE) |> 
  dplyr::select(-month)

# Grouped by location
pendant |> 
  filter(str_detect(season, "wi")) |>
  group_by(location) |> 
  summarise_if(is.numeric, mean, na.rm = TRUE) |> 
  dplyr::select(-month)

pendant |> 
  filter(str_detect(season, "su")) |> 
  group_by(location) |> 
  summarise_if(is.numeric, mean, na.rm = TRUE) |> 
  dplyr::select(-month)

# Urchin species
urchin_orig |> pull(species) |> unique()

urchin_orig |> filter(str_detect(species, "tpil"))

## Transect depths 
df_depth |> 
  mutate(datetime = floor_date(datetime, "month")) |> 
  group_by(location, trans) |> 
  summarise(mean_depth = mean(depth), sd = sd(depth)) 

# Average depth of transects at both locations
df_depth |> 
  mutate(datetime = floor_date(datetime, "month")) |> 
  group_by(trans) |> 
  summarise(mean_depth = mean(depth), sd = sd(depth)) 

# Urchins per square meter
u_sqm = df_urchin |> 
  dplyr::select(datetime, location, trans, species) |> 
  mutate(species = factor(species, levels = c("dsav", "dset", "hcra"))) |> 
  nest(data = species) |> 
  mutate(sp = purrr::map(data, function(x) {
    table(x) |> as_tibble() |> 
      rename(species = species, count = n)
  })) |> unnest(sp) |> dplyr::select(-data) |> 
  arrange(datetime) |> group_by(datetime, location, trans) |> 
  mutate(total = sum(count), 
         percent = count/total*100, 
         sqm = count/20,
         total_per_sqm = total/20) |> 
  ungroup()

# Count of the rugosity measurements
df_rugosity |> 
  group_by(location, trans) |> 
  count()

####################### BRMS global modeling parameters ########################
CHAINS  = 4
CORES   = CHAINS
SEED    = 2021
ITER    = 4000
CONTROL = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.001)

## Model: Wave height ----------------------------------------------------------
waves = df_wave |> 
  mutate(month = month(datetime),
         month = (month - 1)/12,
         year  = year(datetime),
         time = year + month) |> 
  filter(datetime > "2020-08-01")

m_wave_null = bf(mean_height ~ 1) + Gamma(link = "log")

m_wave = bf(mean_height ~ s(time, bs = "cr", k = 6) + location) + Gamma(link = "log")

get_prior(m_wave, data = waves)

waves |> 
  pull(mean_height) |> range()

prior_wn = 
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("student_t(3, 0, 1)", class = "sigma")

prior_w = 
  set_prior("normal(0, 0.5)",     class = "b") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("student_t(3, 0, 1)", class = "sds") +
  set_prior("gamma(0.01, 0.01)",  class = "shape")

m_wave_null_out = 
  brm(m_wave_null,
      data      = waves,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_wn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "transect_wave_null" # Creates a .rds file from the model fit 
  )

m_wave_out = 
  brm(m_wave,
      data      = waves,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_w,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "transect_wave" # Creates a .rds file from the model fit
      )

pp_check(m_wave_out, type = "hist", ndraws = 100)
pp_check(m_wave_out, type = "ecdf_overlay", ndraws = 100)
wave_null = loo(m_wave_null_out)
wave_mod  = loo(m_wave_out)
loo_compare(wave_mod, wave_null)
mcmc_parcoord(m_wave_out)

## Model: Light (PPFD) ----------------------------------------------------------------
# Since the light and temp data has missing values from Jul-Dec 2021,I expanded
# the original data to produce NA's and used the model fit to impute missing data
light = df_light |> expand(datetime, location, trans) |> 
  left_join(df_light, by = c("datetime", "location", "trans")) |> 
  mutate(month = month(datetime),
         month = (month - 1)/12,
         year  = year(datetime),
         time = year + month) |> 
  mutate(datetime = as.Date(datetime)) |> 
  filter(datetime > "2020-08-01")

m_light_null = 
  bf(mean_ppfd ~ 1) + gaussian(link = "log")

m_light = 
  bf(mean_ppfd ~ s(time, bs = "cr", k = 5) + 
               interaction(location, trans) + location + trans) + gaussian(link = "log")

get_prior(m_light, data = light)

prior_ln = 
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("student_t(3, 0, 1)", class = "sigma")

prior_l = 
  set_prior("normal(0, 0.5)", class = "b") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("student_t(3, 0, 1)", class = "sds") +
  set_prior("student_t(3, 0, 1)", class = "sigma")

m_light_null_out = 
  brm(m_light_null,
      data      = light,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_ln,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "transect_ppfd_null" # Creates a .rds file from the model fit
  )

m_light_out = 
  brm(m_light,
      data      = light,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_l,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "transect_ppfd" # Creates a .rds file from the model fit
      )

pp_check(m_light_out, ndraws = 100)
pp_check(m_light_null_out, ndraws = 100)
pp_check(m_light_out, type = "ecdf_overlay", ndraws = 100)
mcmc_rank_hist(m_light_out)
light_null = loo(m_light_null_out)
light_mod  = loo(m_light_out)
loo_compare(light_mod, light_null)

## Model: Temperature ----------------------------------------------------------
# Since the light and temp data has missing values from Jul-Dec 2021, 
# I expanded the original data to produce NA's and used the model fit to impute missing data
temperature = df_light |> 
  expand(datetime, location, trans) |> 
  left_join(df_light, by = c("datetime", "location", "trans")) |> 
  mutate(month = month(datetime),
         month = (month - 1)/12,
         year  = year(datetime),
         time = year + month) |> 
  mutate(datetime = as.Date(datetime)) |> 
  filter(datetime > "2020-08-01")

m_temp_null = 
  bf(mean_temp ~ 1) + gaussian()

m_temp = 
  bf(mean_temp ~ 
       s(time, bs = "cr", k = 12) +
       interaction(location, trans) + location + trans) + gaussian()

get_prior(m_temp, data = temperature)

prior_tn = 
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("student_t(3, 0, 1)", class = "sigma")

prior_t = 
  set_prior("normal(0, 0.5)",     class = "b") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("student_t(3, 0, 1)", class = "sds") +
  set_prior("student_t(3, 0, 1)", class = "sigma")

m_temp_null_out = 
  brm(m_temp_null,
      data      = temperature,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_tn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "transect_temp_null" # Creates a .rds file from the model fit
  )

m_temp_out = 
  brm(m_temp,
      data      = temperature,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_t,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "transect_temp" # Creates a .rds file from the model fit
      )

pp_check(m_temp_out, type = "hist", ndraws = 100)
pp_check(m_temp_out, type = "ecdf_overlay", ndraws = 100)
mcmc_rank_hist(m_temp_out)
loo(m_temp_null_out)
loo(m_temp_out)
loo(m_temp_out, moment_match = TRUE)
loo_compare(temp_mod, temp_null)

## Model: Benthic rugosity -----------------------------------------------------
df_rugosity |> 
  dplyr::select(survey, location, rugosity)

m_rug_null = bf(rugosity ~ 1) + Gamma(link = "log")

m_rug = bf(rugosity ~ location * trans) + Gamma(link = "log")

get_prior(m_rug, data = df_rugosity)

df_rugosity |> 
  print(n = Inf)

prior_rugn = 
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "shape")

prior_rug = 
  set_prior("normal(0, 0.5)",     class = "b") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "shape")

m_rug_null_out = 
  brm(m_rug_null,
      data      = df_rugosity,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_rugn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "transect_rugosity_null" # Creates a .rds file from the model fit
  )

m_rug_out = 
  brm(m_rug,
      data      = df_rugosity,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_rug,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "transect_rugosity" # Creates a .rds file from the model fit
  )

pp_check(m_rug_out, type = "hist", ndraws = 100)
pp_check(m_rug_out, type = "ecdf_overlay", ndraws = 100)
mcmc_rank_hist(m_rug_out)
rug_mod  = loo(m_rug_out)
rug_null = loo(m_rug_null_out)
loo_compare(rug_mod, rug_null)

## Model: Benthic cover --------------------------------------------------------
# I added a row for 2021-01-01 since there was no field work at that time.
# Missing data will be imputed later using GAM
benth = df_quad |> 
  add_row(datetime = ymd_hms("2021-01-01 00:00:00"),
          location = c("Vegetated", "Isoyake"),
          trans    = c("Deep", "Shallow")) |>
  expand(datetime, location, trans, element) |> 
  arrange(datetime) |> 
  drop_na() |> 
  left_join(df_quad, by = c("datetime", "location", "trans", "element")) |>
  filter(!str_detect(element, "Oth")) |> 
  mutate(element = recode(element, 
                          Coralline = "cor", Turf = "tur", Macroalgae = "mac",
                          Substrate = "sub")) |> 
  mutate(month = month(datetime),
         month = (month - 1)/12,
         year = year(datetime),
         time = year + month) |> 
  mutate(datetime = as.Date(datetime)) |> 
  filter(datetime > "2020-08-01") |> 
  dplyr::select(-c(sum, total, month, year, x)) |> 
  pivot_wider(names_from = element, values_from = rate)

# Each benthic element was analyzed separately since the response variables 
# had different family distributions. 
# For example, some benthic elements had excess zeroes while others did not
m_benth_c_null = bf(cor ~ 1) + zero_inflated_beta(link = "logit")
m_benth_m_null = bf(mac ~ 1) + zero_inflated_beta(link = "logit")
m_benth_s_null = bf(sub ~ 1) + zero_one_inflated_beta(link = "logit")
m_benth_t_null = bf(tur ~ 1) + zero_inflated_beta(link = "logit")

m_benth_c = bf(cor ~ 
                 s(time, bs = "cr", k = 5, by = interaction(location, trans)) +
                 interaction(location, trans) + location + trans,
               zi ~ 
                 s(time, bs = "cr", k = 5, by = interaction(location, trans)) +
                 interaction(location, trans) + location + trans) + zero_inflated_beta(link = "logit")

m_benth_m = bf(mac ~ 
                 s(time, bs = "cr", k = 5, by = interaction(location, trans)) +
                 interaction(location, trans) + location + trans,
               zi ~  
                 s(time, bs = "cr", k = 5, by = interaction(location, trans)) +
                 interaction(location, trans) + location + trans) + zero_inflated_beta(link = "logit")

m_benth_s = bf(sub ~ 
                 s(time, bs = "cr", k = 5, by = interaction(location, trans)) +
                 interaction(location, trans) + location + trans) + zero_one_inflated_beta(link = "logit")

m_benth_t = bf(tur ~ 
                 s(time, bs = "cr", k = 5) +
                 interaction(location, trans) + location + trans,
               zi ~ 
                 s(time, bs = "cr", k = 5) +
                 interaction(location, trans) + location + trans) + zero_inflated_beta(link = "logit")

get_prior(m_benth_s, data = benth)

benth |> print(n = Inf)

prior_bn = 
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "phi") +
  set_prior("beta(1, 1)",         class = "zi")

prior_b1n = 
  set_prior("beta(1, 1)",         class = "coi") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "phi")+
  set_prior("beta(1, 1)",         class = "zoi")


prior_b = 
  set_prior("normal(0, 0.5)",     class = "b") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "phi") +
  set_prior("student_t(3, 0, 1)", class = "sds") +
  set_prior("normal(0, 0.5)",     class = "b",         dpar = "zi") +
  set_prior("logistic(0, 1)",     class = "Intercept", dpar = "zi") +
  set_prior("student_t(3, 0, 1)", class = "sds",       dpar = "zi")

prior_b1 = 
  set_prior("normal(0, 0.5)",     class = "b") +
  set_prior("beta(1, 1)",         class = "coi") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "phi") +
  set_prior("student_t(3, 0, 1)", class = "sds") +
  set_prior("beta(1, 1)",         class = "zoi") 


m_benth_c_null_out = 
  brm(m_benth_c_null,
      data      = benth,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_bn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_coralline_null" # Creates a .rds file from the model fit
  )

m_benth_m_null_out = 
  brm(m_benth_m_null,
      data      = benth,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_bn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_macroalgae_null" # Creates a .rds file from the model fit
  )

m_benth_s_null_out = 
  brm(m_benth_s_null,
      data      = benth,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_b1n,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_substrate_null" # Creates a .rds file from the model fit
  )

m_benth_t_null_out = 
  brm(m_benth_t_null,
      data      = benth,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_bn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_turf_null" # Creates a .rds file from the model fit
  )

m_benth_c_out = 
  brm(m_benth_c,
      data      = benth,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_b,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_coralline" # Creates a .rds file from the model fit
      )

m_benth_m_out = 
  brm(m_benth_m,
      data      = benth,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_b,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_macroalgae" # Creates a .rds file from the model fit
      )

m_benth_s_out = 
  brm(m_benth_s,
      data      = benth,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_b1,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_substrate" # Creates a .rds file from the model fit
      )

m_benth_t_out = 
  brm(m_benth_t,
      data      = benth,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_b,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_turf" # Creates a .rds file from the model fit
      )

pp_check(m_benth_c_out, type = "hist", ndraws = 100)
pp_check(m_benth_m_out, type = "hist", ndraws = 100)
pp_check(m_benth_s_out, type = "hist", ndraws = 100)
pp_check(m_benth_t_out, type = "hist", ndraws = 100)

benth_c_null = loo(m_benth_c_null_out)
benth_m_null = loo(m_benth_m_null_out)
benth_s_null = loo(m_benth_s_null_out)
benth_t_null = loo(m_benth_t_null_out)

benth_c_mod = loo(m_benth_c_out, moment_match = TRUE)
benth_m_mod = loo(m_benth_m_out, moment_match = TRUE)
benth_s_mod = loo(m_benth_s_out, moment_match = TRUE)
benth_t_mod = loo(m_benth_t_out, moment_match = TRUE)

loo(m_benth_c_out)
loo(m_benth_m_out)
loo(m_benth_s_out)
loo(m_benth_t_out)

loo_compare(benth_c_mod, benth_c_null)
loo_compare(benth_m_mod, benth_m_null)
loo_compare(benth_s_mod, benth_s_null)
loo_compare(benth_t_mod, benth_t_null)

## Model: Urchin density -------------------------------------------------------

## URCHIN COUNTS
sqm = df_urchin |> 
  dplyr::select(datetime, location, trans, species) |> 
  mutate(species = factor(species, levels = c("dsav", "dset", "hcra"))) |> 
  nest(data = species) |> 
  mutate(sp = purrr::map(data, function(x) {
    table(x) |> as_tibble() |> 
      rename(species = species, count = n)
  })) |> unnest(sp) |> dplyr::select(-data) |> 
  arrange(datetime) |> group_by(datetime, location, trans) |> 
# Calculate the % and urchins per square meter in the 20 square meter transect
  mutate(total_c = sum(count), 
         percent = count/total_c*100, 
         sqm = count/20) |> 
  ungroup() 

# I added a row for 2021-01-01 due to missing data.
# Missing data will be imputed later
urchin_sqm = sqm |>  
  add_row(datetime = ymd_hms("2021-01-01 00:00:00"),
          location = c("Vegetated", "Isoyake"),
          trans    = c("Deep", "Shallow")) |>
  expand(datetime, location, trans, species) |>
  drop_na() |>
  left_join(sqm, by = c("datetime", "location", "trans", "species")) |>
  mutate(month = month(datetime),
         month = (month - 1)/12,
         year  = year(datetime),
         time = year + month) |>
  dplyr::select(datetime, location, trans, species, count, time) |> 
  pivot_wider(names_from = species, values_from = count)


count_dsav_mod_null = bf(dsav ~ 1) + hurdle_negbinomial(link = "log")
count_dset_mod_null = bf(dset ~ 1) + hurdle_negbinomial(link = "log")
count_hcra_mod_null = bf(hcra ~ 1) + hurdle_negbinomial(link = "log")

count_dsav_mod = 
  bf(dsav ~ s(time, bs = "cr", k = 5) + 
       s(time, by = interaction(location, trans)) + 
       interaction(location, trans) + location + trans,
     hu ~ s(time, bs = "cr", k = 5) + 
       s(time, by = interaction(location, trans)) + 
       interaction(location, trans) + location + trans) + hurdle_negbinomial("log")

count_dset_mod = 
  bf(dset ~ s(time, bs = "cr", k = 5) + 
       s(time, by = interaction(location, trans)) + 
       interaction(location, trans) + location + trans,
     hu ~ s(time, bs = "cr", k = 5) + 
       s(time, by = interaction(location, trans)) + 
       interaction(location, trans) + location + trans) + hurdle_negbinomial("log")

count_hcra_mod = 
  bf(hcra ~ s(time, bs = "cr", k = 5) + 
       s(time, by = interaction(location, trans)) + 
       interaction(location, trans) + location + trans,
   hu ~ s(time, bs = "cr", k = 5) + 
       s(time, by = interaction(location, trans)) + 
       interaction(location, trans) + location + trans) + hurdle_negbinomial("log")

get_prior(count_dsav_mod_null, data = urchin_sqm)

prior_sqmn = 
  set_prior("beta(1, 1)",         class = "hu") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "shape") 

prior_sqm = 
  set_prior("normal(0, 0.5)", class = "b") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("student_t(3, 0, 1)", class = "sds") +
  set_prior("gamma(0.01, 0.01)", class = "shape") +
  set_prior("normal(0, 0.5)", class = "b", dpar = "hu") +
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu") +
  set_prior("student_t(3, 0, 1)", class = "sds", dpar = "hu") 
  
dsav_sqm_null_out = 
  brm(count_dsav_mod_null,
      data      = urchin_sqm,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_sqmn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_dsav_count_null" # Creates a .rds file from the model fit
  )

dset_sqm_null_out = 
  brm(count_dset_mod_null,
      data      = urchin_sqm,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_sqmn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_dset_count_null" # Creates a .rds file from the model fit
  )

hcra_sqm_null_out = 
  brm(count_hcra_mod_null,
      data      = urchin_sqm,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_sqmn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_hcra_count_null" # Creates a .rds file from the model fit
  )


dsav_sqm_out = 
  brm(count_dsav_mod,
      data      = urchin_sqm,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_sqm,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_dsav_count" # Creates a .rds file from the model fit
  )

dset_sqm_out = 
  brm(count_dset_mod,
      data      = urchin_sqm,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_sqm,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_dset_count" # Creates a .rds file from the model fit
  )

hcra_sqm_out = 
  brm(count_hcra_mod,
      data      = urchin_sqm,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_sqm,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_hcra_count" # Creates a .rds file from the model fit
  )

pp_check(dsav_sqm_out, type = "hist", ndraws = 100)
pp_check(dset_sqm_out, type = "hist", ndraws = 100)
pp_check(hcra_sqm_out, type = "hist", ndraws = 100)

sqm_dsav_null = loo(dsav_sqm_null_out, moment_match = TRUE) 
sqm_dset_null = loo(dset_sqm_null_out, moment_match = TRUE) 
sqm_hcra_null = loo(hcra_sqm_null_out, moment_match = TRUE) 

sqm_dsav_mod = loo(dsav_sqm_out, moment_match = TRUE)
sqm_dset_mod = loo(dset_sqm_out, moment_match = TRUE)
sqm_hcra_mod = loo(hcra_sqm_out, moment_match = TRUE)

loo_compare(sqm_dsav_mod, sqm_dsav_null)
loo_compare(sqm_dset_mod, sqm_dset_null)
loo_compare(sqm_hcra_mod, sqm_hcra_null)


## Model: Urchin biomass -------------------------------------------------------
urchin_bio = df_urchin |> 
  mutate(species = factor(species, levels = c("dsav", "dset", "hcra"))) |>
  group_by(datetime, location, trans, species) |> 
  summarise(weight = mean(weight)) |>
# Calculate urchin biomass per square meter within 20sqm 
  mutate(total_b = sum(weight),
         total_b_sqm = total_b/20,
         bio_sqm = weight/20) 

urchin_full1 = left_join(new_data_count, urchin_bio, by = c("datetime", "location", "trans", "species")) |>  
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) |>
  dplyr::select(-count, -time) |> 
  mutate(datetime = as.Date(datetime)) |> 
  filter(!str_detect(datetime, "2021-01-01")) |> 
  mutate(datetime = as_datetime(datetime)) |>
  right_join(new_data_count, by = c("datetime", "location", "trans", "species")) |>
  arrange(datetime) |> 
  dplyr::select(datetime, location, trans, species, bio_sqm, time) |> 
  pivot_wider(names_from = species, values_from = bio_sqm) |> 
  mutate(month = month(datetime)) 

bio_dsav_mod_null = bf(dsav ~ 1) + hurdle_gamma(link = "log")
bio_dset_mod_null = bf(dset ~ 1) + hurdle_gamma(link = "log")
bio_hcra_mod_null = bf(hcra ~ 1) + hurdle_gamma(link = "log")

bio_dsav_mod = 
  bf(dsav ~ s(time, bs = "cr", k = 5) + 
       interaction(location, trans) + location + trans,
     hu ~ s(time, bs = "cr", k = 5) + 
       interaction(location, trans) + location + trans) + hurdle_gamma(link = "log")

bio_dset_mod = 
  bf(dset ~ s(time, bs = "cr", k = 5) + 
       interaction(location, trans) + location + trans,
     hu ~ s(time, bs = "cr", k = 5) + 
       interaction(location, trans) + location + trans) + hurdle_gamma(link = "log")

bio_hcra_mod = 
  bf(hcra ~ s(time, bs = "cr", k = 5) + 
       interaction(location, trans) + location + trans,
     hu ~ s(time, bs = "cr", k = 5) + 
       interaction(location, trans) + location + trans) + hurdle_gamma(link = "log")

get_prior(bio_dsav_mod, data = urchin_full1)

prior_bn = 
  set_prior("beta(1, 1)",         class = "hu") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.1, 0.1)",    class = "shape") 

prior_b = 
  set_prior("normal(0, 0.5)",     class = "b") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("student_t(3, 0, 1)", class = "sds") +
  set_prior("gamma(0.1, 0.1)",    class = "shape") +
  set_prior("normal(0, 0.5)",     class = "b",         dpar = "hu") +
  set_prior("logistic(0, 1)",     class = "Intercept", dpar = "hu") +
  set_prior("student_t(3, 0, 1)", class = "sds",       dpar = "hu") 

dsav_bio_null_out = 
  brm(bio_dsav_mod_null,
      data      = urchin_full1,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_bn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_dsav_biomass_null" # Creates a .rds file from the model fit
  )

dset_bio_null_out = 
  brm(bio_dset_mod_null,
      data      = urchin_full1,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_bn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_dset_biomass_null" # Creates a .rds file from the model fit
  )

hcra_bio_null_out = 
  brm(bio_hcra_mod_null,
      data      = urchin_full1,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_bn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_hcra_biomass_null" # Creates a .rds file from the model fit
  )

dsav_bio_out = 
  brm(bio_dsav_mod,
      data      = urchin_full1,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_b,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_dsav_biomass" # Creates a .rds file from the model fit
      )

dset_bio_out = 
  brm(bio_dset_mod,
      data      = urchin_full1,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_b,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_dset_biomass" # Creates a .rds file from the model fit
      )

hcra_bio_out = 
  brm(bio_hcra_mod,
      data      = urchin_full1,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_b,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "benthic_hcra_biomass" # Creates a .rds file from the model fit
      )

pp_check(dsav_bio_out, type = "hist", ndraws = 100)
pp_check(dset_bio_out, type = "hist", ndraws = 100)
pp_check(hcra_bio_out, type = "hist", ndraws = 100)

bio_dsav_null = loo(dsav_bio_null_out)
bio_dset_null = loo(dset_bio_null_out)
bio_hcra_null = loo(hcra_bio_null_out)

bio_dsav_mod = loo(dsav_bio_out, moment_match = TRUE)
bio_dset_mod = loo(dset_bio_out, moment_match = TRUE)
bio_hcra_mod = loo(hcra_bio_out, moment_match = TRUE)

loo_compare(bio_dsav_mod, bio_dsav_null)
loo_compare(bio_dset_mod, bio_dset_null)
loo_compare(bio_hcra_mod, bio_hcra_null)

## Model: Urchin microhabitat --------------------------------------------------
tmp = df_urchin |> 
  dplyr::select(datetime, location, trans, species, size, habitat) |> 
  mutate(location = factor(location, levels = c("Isoyake", "Vegetated")),
         trans = factor(trans, levels = c("Deep", "Shallow")),
         species = factor(species, levels = c("dsav", "dset", "hcra")),
         habitat = factor(habitat, 
                          levels = c("burrow", "crevice", "free"),
                          labels = c("bu", "cr", "fr"))) 

tmp_percent = tmp |> 
  group_by(datetime, location, trans) |>
  count(species, habitat) |> 
  ungroup() |> 
  group_by(datetime, location, trans) |> 
  mutate(total = sum(n)) |> 
  mutate(percent = (n/total))

tmp_size = tmp |> 
  group_by(datetime, location, trans, species, habitat) |> 
  summarise(mean_size = mean(size), sd_size = sd(size)) |> 
  mutate(categ = ifelse(mean_size <= 3.9, "s", 
                 ifelse(mean_size > 3.9 & mean_size <= 5.9, "m",
                 ifelse(mean_size > 5.9, "l", "none"))))

tmp_full = tmp |> 
  expand(datetime, location, trans, species, habitat) |> 
  left_join(tmp_size,    by = c("datetime", "location", "trans", "species", "habitat")) |> 
  left_join(tmp_percent, by = c("datetime", "location", "trans", "species", "habitat")) |>
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) |>
  drop_na() |>
  select(-c(sd_size, mean_size)) |> 
  mutate(datetime = as.Date(datetime),
         month = month(datetime),
         month = (month - 1)/12,
         year  = year(datetime),
         time  = year + month)


uhab_null = bf(percent ~ 1) + Beta(link = "logit")

uhab = bf(percent ~ s(time, bs = "cr", k = 15) + 
            s(time, by = interaction(trans, location))    +
            s(time, by = interaction(trans, species))     +
            s(time, by = interaction(species, location))  +
            s(time, by = interaction(species, habitat))   +
            s(time, by = interaction(habitat, location))  +
            s(time, by = interaction(trans, habitat))     +
            s(time, by = interaction(categ, habitat))     +
            s(time, by = interaction(species, categ))     + 
            trans * location * species * habitat * categ) + Beta(link = "logit")

get_prior(uhab, data = tmp_full)

prior_habn = 
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "phi")

prior_hab = 
  set_prior("normal(0, 0.5)",     class = "b") +
  set_prior("student_t(3, 0, 1)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "phi") +
  set_prior("student_t(3, 0, 1)", class = "sds")

uhab_null_out = 
  brm(uhab_null,
      data       = tmp_full,
      chains     = CHAINS,
      cores      = 4,
      seed       = SEED,
      iter       = ITER,
      prior      = prior_habn,
      control    = CONTROL,
      save_pars  = save_pars(all = TRUE),
      # file_refit = "always", 
      file       = "microhabitat_distribution_null" # Creates a .rds file from the model fit
  )

uhab_mod_out = 
  brm(uhab,
      data       = tmp_full,
      chains     = CHAINS,
      cores      = 20,
      seed       = SEED,
      iter       = ITER,
      prior      = prior_hab,
      control    = CONTROL,
      save_pars  = save_pars(all = TRUE),
      # file_refit = "always", 
      file       = "microhabitat_distribution" # Creates a .rds file from the model fit
  )

pp_check(uhab_mod_out, type = "hist", ndraws = 100)
pp_check(uhab_mod_out, type = "ecdf_overlay", ndraws = 100)
uhab_null = loo(uhab_null_out)
uhab_mod = loo(uhab_mod_out) 
loo(uhab_mod_out) 
loo_compare(uhab_mod, uhab_null)
loo_compare(uhab_mod, uhab_null) |> as_tibble()

## Model: Urchin linear Displacement GLM ---------------------------------------
tag_linear = df_tag |> 
  select(survey, location, trial, species, distance) |>
  mutate(survey = factor(survey, levels = c("start", "mid", "end"), 
                         labels = c("0-hr", "12-hr", "24-hr"))) |> 
  drop_na()

m_disp_null = 
  bf(distance ~ 1) + hurdle_gamma(link = "log")

m_disp_mod = 
  bf(distance ~ 
          species * location * survey + (1|trial),
     hu ~ species * location * survey) + hurdle_gamma(link = "log")

get_prior(m_disp_mod, data = tag_linear)

prior_dispn = 
  set_prior("beta(1, 1)",         class = "hu") +
  set_prior("student_t(3, 0, 2)", class = "Intercept") +
  set_prior("gamma(0.01, 0.01)",  class = "shape") 

prior_disp = 
  set_prior("normal(0, 1)",       class = "b") +
  set_prior("student_t(3, 0, 2)", class = "Intercept") +
  set_prior("student_t(3, 0, 2)", class = "sd") +
  set_prior("gamma(0.01, 0.01)",  class = "shape") +
  set_prior("normal(0, 1)",       class = "b",         dpar = "hu") +
  set_prior("logistic(0, 1)",     class = "Intercept", dpar = "hu") 

disp_null_out = 
  brm(m_disp_null,
      data      = tag_linear,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_dispn,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "tag_distance_null" # Creates a .rds file from the model fit
      )

disp_mod_out = 
  brm(m_disp_mod,
      data      = tag_linear,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_disp,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "tag_distance" # Creates a .rds file from the model fit
      )

pp_check(disp_mod_out, type = "hist", binwidth = 2, ndraws = 3) 
pp_check(disp_mod_out, type = "ecdf_overlay", ndraws = 100)
loo(disp_mod_out)
disp_null = loo(disp_null_out) 
disp_mod  = loo(disp_mod_out)  
loo_compare(disp_mod, disp_null) 

plot(disp_mod_out)
mcmc_rank_hist(disp_mod_out, pars = vars("b_Intercept", "b_specieshcra"))
plot(conditional_effects(disp_mod_out))
conditional_effects(disp_mod_out, dpar = "hu")


## Model: Urchin Group size GLM ------------------------------------------------
tag_groupsize = df_tag |> 
  select(survey, location, trial, species, individuals) |>
  filter(str_detect(species, "dset")) |>
  mutate(survey = factor(survey, levels = c("start", "mid", "end"), 
                         labels = c("0-hr", "12-hr", "24-hr"))) |> 
  drop_na()

group_null = bf(individuals ~ 1) + negbinomial(link = "log")

group_model = bf(individuals ~ survey * location + (location + survey|trial)) + negbinomial(link = "log")

get_prior(group_model, data = tag_groupsize)

prior_group0 = 
  set_prior("student_t(3, 0, 2)", class = "Intercept") +
  set_prior("student_t(3, 0, 2)", class = "shape")

prior_group1 = 
  set_prior("student_t(3, 0, 2)", class = "b")         +
  set_prior("student_t(3, 0, 2)", class = "Intercept") +
  set_prior("student_t(3, 0, 2)", class = "sd")        +
  set_prior("student_t(3, 0, 2)", class = "shape")

group_null_out = 
  brm(group_null,
      data      = tag_groupsize,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_group0,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "tag_group_null" # Creates a .rds file from the model fit
  )

group_mod_out = 
  brm(group_model,
      data      = tag_groupsize,
      chains    = CHAINS,
      cores     = CORES,
      seed      = SEED,
      iter      = ITER,
      prior     = prior_group1,
      control   = CONTROL,
      save_pars = save_pars(all = TRUE),
      # file_refit = "always", 
      file      = "tag_group" # Creates a .rds file from the model fit
      )

pp_check(group_mod_out,  type = "hist", ndraws = 5)
pp_check(group_mod_out,  type = "ecdf_overlay", ndraws = 100)
group_null = loo(group_null_out)
group_mod  = loo(group_mod_out)
loo_compare(group_mod, group_null)
plot(conditional_effects(group_mod_out))

