# Creating figures/plots for the manuscript:
# "The behavior of sympatric urchin species across and ecosystem state gradient"

# Load packages for data reading, analysis, and plotting
library(tidyverse)  # For general data wrangling
library(patchwork)  # For easy positioning of multiple plots
library(lubridate)  # For easy manipulation of dates
library(lemon)      # For added ggplot functionality
library(ggpubr)     # For added ggplot functionality
library(scales)     # For added ggplot functionality
library(tidybayes)  # For extracting tidy data from Bayesian model fits
library(bayesplot)  # For added plotting functionality of Bayesian models
library(readxl)     # For reading excel files
library(gnnlab)     # For easy reading of HOBO logger data
library(showtext)   # For specifying text/graphics

# Map making packages
library(magick)     # For general image processing
library(pdftools)   # For 
library(ggspatial)  # For working with spatial data in ggplot
library(sf)         # Functions for working with spatial feature data
library(ggrepel)    # Functions for adding distance between text/labels in ggplot
library(rgdal)      # Functions for working with geospatial data

# Set global parameters for seed, font, text and system locale
set.seed(2021)
font_add_google("Noto Sans","notosans")
theme_pubr(base_size = 10, base_family = "notosans") |> theme_set()
showtext_auto()
Sys.setlocale("LC_TIME", "en_US.UTF-8")
################################ LOAD DATA #####################################
df_urchin   = read_csv("urchin_data.csv")         # Urchin size-weight model fit data
df_quad     = read_csv("transect_quad_data.csv")  # Benthic quadrat data. Excludes "Others"
df_depth    = read_csv("depth_data.csv")          # Depth logger data
df_l_t      = read_csv("temp_ppfd.csv")           # Mean monthly temperature and light (ppfd) data
df_wave     = read_csv("wave_data.csv")           # Wave logger data
df_tag      = read_csv("tag_urchin.csv")          # Mark-recapture experiment data
temp_mean   = read_csv("winter_temp.csv")         # 2018-2021 Vegetated site winter temperatures
df_rugosity = read_csv("raw_rugosity.csv")        # Benthic rugosity data

################################ FIGURES #######################################
# Figure 1(Study site map) -----------------------------------------------------

# Load shapefiles
shpfile  = "N03-20210101_42_GML/N03-21_42_210101.shp"
shpfile2 = "W05-07_42_GML/W05-07_42-g_RiverNode.shp"
shpfile3 = "W05-07_42_GML/W05-07_42-g_Stream.shp"
shpfile4 = "L03-b-u-09_4929-jgd_GML/L03-b-u-09_4929.shp"

# Inset map and feature coordinates 
# Set the transect coordinates
df = read_csv("transect_coords.csv",
              locale = locale(encoding = "sjis")) 
main = df |> 
  filter(str_detect(location, "Ari")) |> 
  mutate(label = " ")


inset = df |>
  filter(!str_detect(location, "Ari")) |> 
  mutate(label = paste(location, trans, sep = ": ")) |> 
  mutate(x_nudge = case_when(label == "Vegetated: Deep" ~ 0.002,
                             label == "Vegetated: Shallow" ~ 0.002,
                             label == "Isoyake: Deep" ~ 0.001,
                             label == "Isoyake: Shallow" ~ 0.002),
         y_nudge = case_when(label == "Vegetated: Deep" ~ 0.0008,
                             label == "Vegetated: Shallow" ~ -0.00005,
                             label == "Isoyake: Deep" ~ 0.0005,
                             label == "Isoyake: Shallow" ~ -0.0004))

# Read shapefiles into objects
kamigoto = st_read(shpfile,  as_tibble = TRUE)
crs_original = st_crs(kamigoto)

water    = st_read(shpfile2, as_tibble = TRUE, crs = crs_original)
stream   = st_read(shpfile3, as_tibble = TRUE, crs = crs_original)
land     = st_read(shpfile4, as_tibble = TRUE, options = "ENCODING=CP932")

# Boundaries for map
lonlimits_map = c(128.557478, 129.291732)
latlimits_map = c(32.517132, 33.336734)

# Boundaries for the map inset
lonlimits = c(129.117251, 129.123130)
latlimits = c(32.990712, 32.983944)

box = c(xmin = lonlimits[1], ymin = latlimits[1],
        xmax = lonlimits[2], ymax = latlimits[2])
water = st_crop(water, box)
stream= st_crop(stream, box)
land  = st_transform(land, crs = crs_original)
land  = st_crop(land, box)
trans_main = st_as_sf(main, coords = c("long", "lat"), crs = crs_original) 
trans_inset= st_as_sf(inset, coords = c("long", "lat"), crs = crs_original) 


ggplot() + 
  geom_sf(data = kamigoto, size = 0.5, alpha = 1, color = "black") +
  geom_text_repel(aes(long, lat, label = label),
                  size = 10,
                  data = inset,
                  color = "black",
                  point.padding = 0.2,
                  nudge_x = inset$x_nudge,
                  nudge_y = inset$y_nudge) +
  geom_sf(aes(color = location), data = trans_inset, size = 0.5, alpha = 0.9) +
  coord_sf(xlim = lonlimits, ylim = latlimits, expand = F) +
  theme_pubr() +
  # Add a scale bar to the map
  annotation_scale(height = unit(0.25, "cm"), text_cex = 2) + 
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

H = 100
W = 100

ggsave(filename = "arikawa_inset.png", 
       width = W, height = H, units = "mm")

# Background map coordinates (Map of the Goto Islands)
fname = "japan_ver81.shp"
jp = readOGR(fname)
jp = fortify(jp) 

# Goto map coordinates/boundaries
longA = 128.557478
longB = 129.291732
latA  = 32.517132
latB  = 33.336734

Fontsize = 18
head(jp)
ggplot() + 
  geom_polygon(aes(long, lat, group = group),
               data = jp,
               fill = "grey90",
               colour = "Black") + 
  coord_map(xlim = c(longA, longB), ylim = c(latA, latB)) +
  geom_text_repel(aes(long, lat, label = label),
                  data = main,
                  color = "black",
                  point.padding = 0.2,
                  nudge_x = 1,
                  nudge_y = -0.09) +
  theme_pubr() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 20),
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_rect(fill = "lightblue2")
  )

H = 100
W = 100

ggsave(filename = "arikawa_map.png", 
       width = W, height = H, units = "mm")

# Combine background and inset 
pic1 = "arikawa_map.png"   |> image_read()
pic2 = "arikawa_inset.png" |> image_read()

img1  = image_resize(pic1,  geometry = "1000x1000!")
img2  = image_resize(pic2,  geometry = "400x400!") 

# template for cropping
crop1_w = image_info(img1) |> pull(width) - 50
crop1_h = image_info(img1) |> pull(height)
crop2_w = image_info(img2) |> pull(width) - 125
crop2_h = image_info(img2) |> pull(height) - 15

# Crop the images using the template
A1 = image_crop(img1, geometry=paste0(crop1_w,"x",crop1_h,"+120+0"),repage=TRUE)
A2 = image_crop(img2, geometry=paste0(crop2_w,"x",crop2_h,"+65+0"),repage=TRUE) |> 
  image_border(color = "black", geometry = "1x1")

# Background
W = image_info(A1) |> pull(width) + 220
H = image_info(A1) |> pull(height)
bg1 = image_blank(W, H, color = "rgba(0, 0, 0, 0") # 100% transparency

image1 = image_composite(A1, bg1, gravity = "center", offset = "+0+0") |> 
  image_border(color = "white", geometry = "1x1")

full_map =
  image_composite(image1, A2, operator = "atop", gravity = "southeast", offset = "+20+110") |> 
  image_border(color = "black", geometry = "1x1")

w = 100
h = 100

image_write(full_map, 
            format = "pdf", path = "./Fig_1.pdf", density = 100)


# Figure 2(Logger data) --------------------------------------------------------

# Wave (load .rds model fits from urchin_ecosystem_analysis.R)
m_wave_null_out = read_rds("transect_wave_null.rds")
m_wave_out      = read_rds("transect_wave.rds")

waves = df_wave |> 
  mutate(month = month(datetime),
         month = month/12,
         year  = year(datetime),
         x = month*12,
         time = year + month) |> 
  filter(datetime > "2020-08-01")

wave_out = m_wave_out$data |> add_epred_draws(m_wave_out, seed = 2022) |> 
  group_by(time, location) |> 
  mean_hdci(.epred) |> 
  left_join(waves, by = c("location", "time"))

wave1 = 
  wave_out |> 
  ggplot() + 
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = location), alpha = 0.3)+
  geom_line(aes(x = datetime, y = .epred, color = location, linetype = location), size = 0.5) +
  geom_point(aes(x = datetime, y = mean_height, color = location), size = 0.5) +
  scale_x_date("Month", minor_breaks = "month", labels = label_date_short()) +
  scale_y_continuous("Wave height (m)") +
  scale_fill_viridis_d( end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  theme_pubr() +
  ggtitle("(A)") +
  theme(
    legend.title      = element_blank(),
    legend.background = element_blank(),
    legend.position   = c(0.78,0.92),
    strip.text        = element_text(color = "white"),
    strip.background  = element_rect(fill = "black"),
    axis.text         = element_text(size = 10),
    axis.text.x       = element_blank(),
    legend.key.size   = unit(0.80, "lines"),
    axis.title.x      = element_blank(),
  ) 

# Light (load .rds model fits from urchin_ecosystem_analysis.R)
m_light_null_out = read_rds("transect_ppfd_null.rds")
m_light_out      = read_rds("transect_ppfd.rds")

light = df_l_t |> expand(datetime, location, trans) |> 
  left_join(df_l_t, by = c("datetime", "location", "trans")) |> 
  mutate(month = month(datetime),
         month = month/12,
         year  = year(datetime),
         x     = month*12,
         time = year + month) |> 
  mutate(datetime = as.Date(datetime)) |> 
  filter(datetime > "2020-08-01")

light_out = add_epred_draws(m_light_out, newdata = light) |> 
  group_by(time, location, trans) |> 
  mean_hdci(.epred) |> 
  left_join(light, by = c("location", "trans", "time"))

ppfd = 
  light_out |> 
  ggplot() + 
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = location), alpha = 0.3)+
  geom_line(aes(x = datetime, y = .epred, color = location, linetype = location), size = 0.5) +
  geom_point(aes(datetime, mean_ppfd, color = location), size = 0.5) +
  scale_y_continuous(expression(paste("Light",~("PPFD, "~mol~m^{-2}~day^{-1})))) +
  scale_x_date("Month", minor_breaks = "month", labels = label_date_short()) +
  facet_rep_grid(rows = vars(trans)) +
  scale_fill_viridis_d(end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  theme_pubr() +
  ggtitle("(B)") +
  theme(
    legend.position  = "none",
    strip.text.x     = element_blank(),
    strip.text       = element_text(color = "white"),
    strip.background = element_rect(fill = "black"),
    axis.text        = element_text(size = 10),
    axis.text.x      = element_blank(),
    axis.title.x     = element_blank()
  ) 

# Temperature (load .rds model fits from urchin_ecosystem_analysis.R)
m_temp_null_out = read_rds("transect_temp_null.rds")
m_temp_out      = read_rds("transect_temp.rds")

temperature = df_l_t |> 
  expand(datetime, location, trans) |> 
  left_join(df_l_t, by = c("datetime", "location", "trans")) |> 
  mutate(month = month(datetime),
         month = month/12,
         year  = year(datetime),
         x     = month*12,
         time  = year + month) |> 
  mutate(datetime = as.Date(datetime)) |> 
  filter(datetime > "2020-08-01")

temp_out = add_epred_draws(m_temp_out, newdata = temperature) |> 
  group_by(time, location, trans) |> 
  mean_hdci(.epred) |> 
  left_join(light, by = c("location", "trans", "time"))

temp = 
  temp_out |> 
  ggplot() + 
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = location), alpha = 0.3)+
  geom_line(aes(x = datetime, y = .epred, color = location, linetype = location), size = 0.5) +
  geom_point(aes(x = datetime, y = mean_temp, color = location), size = 0.5) +
  scale_y_continuous("Temperature (°C)", 
                     limits = c(10, 30)
                     ) +
  scale_x_date("Month", minor_breaks = "month", labels = label_date_short()) +
  facet_rep_grid(rows = vars(trans)) +
  scale_fill_viridis_d(end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  theme_pubr() +
  ggtitle("(C)") +
  theme(
    legend.position  = "none",
    strip.text.x     = element_blank(),
    strip.text       = element_text(color = "white"),
    strip.background = element_rect(fill = "black"),
    panel.spacing.y  = unit(-1, "lines"),
    axis.text        = element_text(size = 10)
  )

w = 100
h = 220

wave1/ppfd/temp

ggsave("Fig_2.pdf", width = w, height = h, units = "mm")

wave_null  = loo(m_wave_null_out)
light_null = loo(m_light_null_out)
temp_null  = loo(m_temp_null_out)
wave_mod   = loo(m_wave_out)
light_mod  = loo(m_light_out)
temp_mod   = loo(m_temp_out, moment_match = TRUE)


# Figure 3(Rugosity) -----------------------------------------------------------
# (load .rds model fits from urchin_ecosystem_analysis.R)
m_rug_null_out = read_rds("transect_rugosity_null.rds")
m_rug_out      = read_rds("transect_rugosity.rds")

rugosity = df_rugosity |> 
  dplyr::select(survey, location, trans, rugosity) |> 
  mutate(trans = recode(trans, deep = "Deep", shallow = "Shallow"))

m_rug_out$data |> 
  add_epred_draws(m_rug_out, seed = 2022) |> 
  group_by(location, trans) |> 
  mean_hdci(.epred) |> 
  mutate(trans = recode(trans, deep = "Deep", shallow = "Shallow")) |> 
  ggplot(aes(location, .epred)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper, color = trans),
                  position = position_dodge(width = 0.5)) +
  geom_point(aes(location, rugosity, color = trans),
             alpha = 0.2,
             data = rugosity,
             position = position_dodge(width = 0.2)) +
  theme_pubr() +
  scale_y_continuous("Rugosity") +
  xlab("Habitat") +
  scale_color_viridis_d(end = 0.7) +
  theme(
    legend.title      = element_blank(),
    legend.background = element_blank(),
    legend.position   = c(0.15, 0.9)
  )

w = 100
h = 100

ggsave(path = "folder_figures", "Fig_3.pdf", width = w, height = h, units = "mm")

rug_mod  = loo(m_rug_out)
rug_null = loo(m_rug_null_out)

# Figure 4(Benthos) ------------------------------------------------------------
# (load .rds model fits from urchin_ecosystem_analysis.R)
m_benth_c_null_out = read_rds("benthic_coralline_null.rds")
m_benth_m_null_out = read_rds("benthic_macroalgae_null.rds")
m_benth_s_null_out = read_rds("benthic_substrate_null.rds")
m_benth_t_null_out = read_rds("benthic_turf_null.rds")

m_benth_c_out = read_rds("benthic_coralline.rds")
m_benth_m_out = read_rds("benthic_macroalgae.rds")
m_benth_s_out = read_rds("benthic_substrate.rds")
m_benth_t_out = read_rds("benthic_turf.rds")

benth = df_quad |> 
  add_row(datetime = ymd_hms("2021-01-01 00:00:00"),
          location = c("Vegetated", "Isoyake"),
          trans    = c("Deep", "Shallow")) |>
  expand(datetime, location, trans, element) |> 
  arrange(datetime) |> 
  drop_na() |> 
  left_join(df_quad, by = c("datetime", "location", "trans", "element")) |> 
  mutate(element = recode(element, 
                          Coralline = "cor", Turf = "tur", Macroalgae = "mac",
                          Substrate = "sub")) |> 
  mutate(month = month(datetime),
         month = month/12,
         year  = year(datetime),
         x     = month*12,
         time  = year + month) |> 
  mutate(datetime = as.Date(datetime)) |> 
  filter(datetime > "2020-08-01") |> 
  dplyr::select(-c(sum, total, month, year, x)) |> 
  pivot_wider(names_from = element, values_from = rate)

new_dat = benth |>
  pivot_longer(cols = cor:tur, names_to = "element", values_to = "rate")

cor_out = add_epred_draws(m_benth_c_out, newdata = benth |> dplyr::select(-c(mac, sub, tur)), seed = 2022) |> 
  group_by(time, location, trans) |> 
  mean_hdci(.epred) |> 
  mutate(element = rep_len("cor", length.out = n()))

mac_out = add_epred_draws(m_benth_m_out, newdata = benth |> dplyr::select(-c(cor, sub, tur)), seed = 2022) |> 
  group_by(time, location, trans) |> 
  mean_hdci(.epred) |>
  mutate(element = rep_len("mac", length.out = n()))

sub_out = add_epred_draws(m_benth_s_out, newdata = benth |> dplyr::select(-c(cor, mac, tur)), seed = 2022) |> 
  group_by(time, location, trans) |> 
  mean_hdci(.epred) |>
  mutate(element = rep_len("sub", length.out = n()))

tur_out = add_epred_draws(m_benth_t_out, newdata = benth |> dplyr::select(-c(cor, mac, sub)), seed = 2022) |> 
  group_by(time, location, trans) |> 
  mean_hdci(.epred) |>
  mutate(element = rep_len("tur", length.out = n()))

benth_out = bind_rows(cor_out, mac_out, sub_out, tur_out) |> 
  left_join(new_dat, by = c("location", "trans", "time", "element")) |> 
  mutate(element = factor(element, 
                          levels = c("mac", "cor", "tur", "sub"),
                          labels = c("Macroalgae", "Coralline", "Turf", "Substrate"))) |> 
  mutate_at(c(".epred", ".lower", ".upper", "rate"), ~ .*100)

benth_out |> print(n = Inf)

benthos =
  benth_out |> 
  ggplot() + 
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = element), alpha = 0.3)+
  geom_line(aes(x = datetime, y = .epred, linetype = element), size = 0.3) +
  geom_point(aes(x = datetime, y = rate, color = element), size = 0.3) +
  guides(linetype = "none") +
  scale_fill_viridis_d(end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  scale_y_continuous("Benthic cover (%)") +
  scale_x_date("Month", date_minor_breaks = "month", labels = label_date_short()) +
  facet_rep_grid(cols = vars(location), rows = vars(trans)) +
  ggtitle("(A)") +
  theme_pubr() +
  theme(
    legend.position   = c(0.08, 0.80),
    legend.key.size   = unit(0.30, "lines"),
    legend.background = element_blank(),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8),
    legend.spacing.x  = unit(0.3, "mm"),
    strip.text        = element_text(color = "white"),
    strip.background  = element_rect(fill = "black"),
    axis.text.x       = element_blank(),
    axis.title.x      = element_blank(),
    axis.text         = element_text(size = 10)
  )

# Urchin density sqm (load .rds model fits from urchin_ecosystem_analysis.R)
dsav_sqm_null_out = read_rds("benthic_dsav_count_null.rds")
dset_sqm_null_out = read_rds("benthic_dset_count_null.rds")
hcra_sqm_null_out = read_rds("benthic_hcra_count_null.rds")

dsav_sqm_out = read_rds("benthic_dsav_count.rds")
dset_sqm_out = read_rds("benthic_dset_count.rds")
hcra_sqm_out = read_rds("benthic_hcra_count.rds")

sqm = df_urchin |> 
  dplyr::select(datetime, location, trans, species) |> 
  mutate(species = factor(species, levels = c("dsav", "dset", "hcra"))) |> 
  nest(data = species) |> 
  mutate(sp = purrr::map(data, function(x) {
    table(x) |> as_tibble() |> 
      rename(species = species, count = n)
  })) |> unnest(sp) |> dplyr::select(-data) |> 
  arrange(datetime) |> group_by(datetime, location, trans) |> 
  mutate(total_c = sum(count), percent = count/total_c*100, sqm = count/20) |> 
  ungroup() |> 
# I added a row for 2021-01-01 since there was no field work at that time.
  add_row(datetime = ymd_hms("2021-01-01 00:00:00"),
          location = c("Vegetated", "Isoyake"),
          trans    = c("Deep", "Shallow"))


# Missing data will be imputed later using GAM
urchin_sqm = sqm |>  
  expand(datetime, location, trans, species) |>
  drop_na() |>
  left_join(sqm, by = c("datetime", "location", "trans", "species")) |>
  mutate(month = month(datetime),
         month = month/12,
         year  = year(datetime),
         x     = month*12,
         time  = year + month) |>
  dplyr::select(datetime, location, trans, species, count, time) |> 
  pivot_wider(names_from = species, values_from = count)

new_data_count = urchin_sqm |> 
  pivot_longer(cols = dsav:hcra, names_to = "species", values_to = "count") 

dsav_count = 
  add_epred_draws(dsav_sqm_out, newdata = urchin_sqm |> dplyr::select(-c(dset, hcra)), seed = 2022) |> 
  group_by(time, location, trans) |> 
  mean_hdci(.epred) |>
  mutate(species = rep_len("dsav", length.out = n()))

dset_count = 
  add_epred_draws(dset_sqm_out, newdata = urchin_sqm |> dplyr::select(-c(dsav, hcra)), seed = 2022) |> 
  group_by(time, location, trans) |> 
  mean_hdci(.epred) |>
  mutate(species = rep_len("dset", length.out = n()))

hcra_count = 
  add_epred_draws(hcra_sqm_out, newdata = urchin_sqm |> dplyr::select(-c(dsav, dset)), seed = 2022) |> 
  group_by(time, location, trans) |> 
  mean_hdci(.epred) |>
  mutate(species = rep_len("hcra", length.out = n()))

urchin_count = bind_rows(dsav_count, dset_count, hcra_count) |> 
  left_join(new_data_count, by = c("location", "trans", "time", "species")) |>
# Convert counts and model fit to urchins/m^2
  mutate_at(c(".epred", ".lower", ".upper", "count"), ~./20) |> 
  mutate(species = recode(species, 
                          dsav = "D. savignyi", 
                          dset = "D. setosum",
                          hcra = "H. crassispina")) |> 
  mutate(datetime = as.Date(datetime))

urchin_count |> print(n = Inf)

urchin_1 =
  urchin_count |> 
  ggplot() + 
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = species), alpha = 0.3) +
  geom_line(aes(x = datetime, y = .epred, linetype = species), size = 0.3) +
  geom_point(aes(x = datetime, y = count, color = species), size = 0.3) +
  guides(linetype = "none") + # This removes the linetype aesthetic to allow italics
  scale_fill_viridis_d(end = 0.8, 
                       option = "B",
                       labels = c(expression(italic("D. savignyi")),
                                  expression(italic("D. setosum")),
                                  expression(italic("H. crassispina")))) +
  scale_color_viridis_d(end = 0.8,
                        option = "B",
                        labels = c(expression(italic("D. savignyi")),
                                   expression(italic("D. setosum")),
                                   expression(italic("H. crassispina")))) +
  scale_y_continuous(expression(paste("Urchin density " (~indiv. ~m^{-2})))) +
  scale_x_date("Month", date_minor_breaks = "month", labels = label_date_short()) +
  facet_rep_grid(cols = vars(location), rows = vars(trans)) +
  theme_pubr() +
  ggtitle("(B)") +
  theme(
    legend.key.size   = unit(0.30, "lines"),
    legend.background = element_blank(),
    legend.title      = element_blank(),
    legend.position   = c(0.095, 0.94),
    legend.text       = element_text(size = 8),
    legend.spacing.x  = unit(0.3, "mm"),
    legend.text.align = 0,
    strip.text.x      = element_blank(),
    strip.text        = element_text(color = "white"),
    strip.background  = element_rect(fill = "black"),
    axis.text.x       = element_blank(),
    axis.text         = element_text(size = 10),
    axis.title.x      = element_blank()
  )

# Urchin biomass sqm (load .rds model fits from urchin_ecosystem_analysis.R)
dsav_bio_null_out = read_rds("benthic_dsav_biomass_null.rds")
dset_bio_null_out = read_rds("benthic_dset_biomass_null.rds")
hcra_bio_null_out = read_rds("benthic_hcra_biomass_null.rds")

dsav_bio_out = read_rds("benthic_dsav_biomass.rds")
dset_bio_out = read_rds("benthic_dset_biomass.rds")
hcra_bio_out = read_rds("benthic_hcra_biomass.rds")

urchin_bio = df_urchin |> 
  mutate(species = factor(species, levels = c("dsav", "dset", "hcra"))) |>
  group_by(datetime, location, trans, species) |> 
  summarise(weight = mean(weight)) |> 
# biomass per square meter within 20sqm 
  mutate(total_b     = sum(weight),
         total_b_sqm = total_b/20,
         bio_sqm     = weight/20) 

urchin_full1 = left_join(new_data_count, urchin_bio, by = c("datetime", "location", "trans", "species")) |>  
  dplyr::select(-count, -time) |> 
  mutate(datetime = as.Date(datetime)) |> 
  right_join(new_data_count, by = c("datetime", "location", "trans", "species")) |>
  arrange(datetime) |> 
  dplyr::select(datetime, location, trans, species, bio_sqm, time) |> 
  pivot_wider(names_from = species, values_from = bio_sqm) |> 
  mutate(month = month(datetime))

new_data_biomass = 
  urchin_full1 |> pivot_longer(cols = dsav:hcra, names_to = "species" , values_to = "bio_sqm")

dsav_biomass = 
  add_epred_draws(dsav_bio_out, newdata = urchin_full1 |> dplyr::select(-c(dset, hcra)), seed = 2022) |>  
  group_by(location, trans, time)  |>  
  mean_hdci(.epred) |> 
  mutate(species = rep_len("dsav", length.out = n())) 

dset_biomass = 
  add_epred_draws(dset_bio_out, newdata = urchin_full1 |> dplyr::select(-c(dsav, hcra)), seed = 2022) |>  
  group_by(location, trans, time) |> 
  mean_hdci(.epred)  |>  
  mutate(species = rep_len("dset", length.out = n())) 

hcra_biomass = 
  add_epred_draws(hcra_bio_out, newdata = urchin_full1 |> dplyr::select(-c(dsav, dset)), seed = 2022) |>  
  group_by(location, trans, time)  |>  
  mean_hdci(.epred)  |>  
  mutate(species = rep_len("hcra", length.out = n())) 

urchin_biomass = bind_rows(dsav_biomass, dset_biomass, hcra_biomass) |> 
  left_join(new_data_biomass, by = c("location", "trans", "species", "time")) |>
  mutate(species = recode(species,
                          dsav = "D. savignyi",
                          dset = "D. setosum",
                          hcra = "H. crassispina")) |> 
  mutate(datetime = as.Date(datetime))

urchin_2 =
  urchin_biomass |> 
  ggplot() +
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = species), alpha = 0.25) +
  geom_line(aes(datetime, .epred, linetype = species), size = 0.3) +
  geom_point(aes(datetime, bio_sqm, color = species), size = 0.5) +
  guides(linetype = "none") +
  facet_rep_grid(cols = vars(location), rows = vars(trans)) +
  scale_x_date("Month", date_minor_breaks = "month", labels = label_date_short()) +
  theme_pubr() +
  scale_y_continuous(expression(paste("Urchin biomass "("grams "~m^{-2})))) +
  ylab("Urchin biomass ") +
  scale_fill_viridis_d (end = 0.8, option = "B") +
  scale_color_viridis_d(end = 0.8, option = "B") +
  ggtitle("(C)") +
  theme(
    legend.position  = "none",
    strip.text.x     = element_blank(),
    strip.text       = element_text(color = "white"),
    strip.background = element_rect(fill = "black"),
    panel.spacing.y  = unit(-1, "lines"),
    axis.text        = element_text(size = 10)
  )

w = 150
h = 220

benthos/urchin_1/urchin_2

ggsave("Fig_4.pdf", width = w, height = h, units = "mm")

benth_c_null = loo(m_benth_c_null_out)
benth_m_null = loo(m_benth_m_null_out)
benth_s_null = loo(m_benth_s_null_out)
benth_t_null = loo(m_benth_t_null_out)
benth_c_mod  = loo(m_benth_c_out, moment_match = TRUE)
benth_m_mod  = loo(m_benth_m_out, moment_match = TRUE)
benth_s_mod  = loo(m_benth_s_out, moment_match = TRUE)
benth_t_mod  = loo(m_benth_t_out, moment_match = TRUE)

sqm_dsav_null = loo(dsav_sqm_null_out, moment_match = TRUE) 
sqm_dset_null = loo(dset_sqm_null_out, moment_match = TRUE) 
sqm_hcra_null = loo(hcra_sqm_null_out, moment_match = TRUE) 
sqm_dsav_mod  = loo(dsav_sqm_out, moment_match = TRUE)
sqm_dset_mod  = loo(dset_sqm_out, moment_match = TRUE)
sqm_hcra_mod  = loo(hcra_sqm_out, moment_match = TRUE)

bio_dsav_null = loo(dsav_bio_null_out)
bio_dset_null = loo(dset_bio_null_out)
bio_hcra_null = loo(hcra_bio_null_out)
bio_dsav_mod  = loo(dsav_bio_out, moment_match = TRUE)
bio_dset_mod  = loo(dset_bio_out, moment_match = TRUE)
bio_hcra_mod  = loo(hcra_bio_out, moment_match = TRUE)

# Figure 5(Microhabitat GAM) ---------------------------------------------------
# (load .rds model fits from urchin_ecosystem_analysis.R)
uhab_mod_null_out = read_rds("microhabitat_distribution_null.rds")
uhab_mod_out      = read_rds("microhabitat_distribution.rds")

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
  dplyr::select(-c(sd_size, mean_size)) |> 
  mutate(datetime = as.Date(datetime),
         month = month(datetime),
         month = month/12,
         year  = year(datetime),
         x = month*12,
         time = year + month) |> 
  dplyr::select(datetime, location, trans, species, habitat, categ, time, percent)

uhab_tmp = uhab_mod_out$data |> add_epred_draws(uhab_mod_out, seed = 2022) |> 
  group_by(time, location, trans, species, habitat, categ) |> 
  mean_hdci(.epred) |> 
  left_join(tmp_full, by = c("time", "location", "trans", "species", "habitat", "categ")) |> 
  mutate_at(c(".epred", ".lower", ".upper", "percent"), ~.*100) |>
  mutate(habitat = recode(habitat, bu = "Pit", cr = "Crevice", fr = "Free-living")) |> 
  mutate(species = recode(species, 
                          dsav = "D. savignyi", 
                          dset = "D. setosum",
                          hcra = "H. crassispina")) |> 
  mutate(categ = factor(categ,
                        levels = c("s", "m", "l"),
                        labels = c("Small (1-3 cm)", "Medium (4-5 cm)", "Large (> 5 cm)")))

uhab1 = uhab_tmp |> 
  filter(str_detect(location, "Iso")) |> 
  filter(str_detect(trans, "D")) |> 
  ggplot() +
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = categ), alpha = 0.3) +
  geom_line(aes(x = datetime, .epred, linetype = categ), lineend = "round") +
  geom_point(aes(datetime, percent, color = categ)) +
  facet_rep_grid(cols = vars(species), rows = vars(habitat)) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_x_date("Month", date_minor_breaks = "month", labels = label_date_short()) +
  scale_y_continuous("Occurrence rate", limits = c(0, 100)) +
  ggtitle("(A) Isoyake: Deep") +
  theme_pubr() +
  scale_fill_viridis_d (end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  theme(
    plot.title        = element_text(size = rel(1.6)),
    legend.title      = element_blank(),
    legend.text       = element_text(size = rel(1)),
    legend.position   = c(0.12, 0.945),
    legend.key.size   = unit(0.8, "lines"),
    legend.background = element_blank(),
    strip.background  = element_rect(fill = "black"),
    strip.text.x      = element_text(face = "italic"),
    strip.text        = element_text(color = "white", size = rel(1.2)),
    panel.spacing.x   = unit(-1.5, "lines"),
    axis.title.x      = element_blank(),
    axis.text.x       = element_blank(),
    axis.text.y       = element_text(size = rel(1)),
    axis.title.y      = element_text(size = rel(1.5))
  )

uhab2 = uhab_tmp |> 
  filter(str_detect(location, "Iso")) |> 
  filter(str_detect(trans, "S")) |> 
  ggplot() +
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = categ), alpha = 0.3) +
  geom_line(aes(x = datetime, .epred, linetype = categ), lineend = "round") +
  geom_point(aes(datetime, percent, color = categ)) +
  facet_rep_grid(cols = vars(species), rows = vars(habitat)) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_x_date("Month", date_minor_breaks = "month", labels = label_date_short()) +
  scale_y_continuous("Occurrence rate", limits = c(0, 100)) +
  ggtitle("(B) Isoyake: Shallow") +
  theme_pubr() +
  scale_fill_viridis_d (end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  theme(
    plot.title       = element_text(size = rel(1.6)),
    legend.position  = "none",
    strip.background = element_rect(fill = "black"),
    strip.text.x     = element_blank(),
    strip.text       = element_text(color = "white", size = rel(1.2)),
    axis.text.y      = element_text(size = rel(1)),
    axis.title       = element_text(size = rel(1.5)),
    panel.spacing.x  = unit(-1.5, "lines"),
    panel.spacing.y  = unit(-1.5, "lines")
  )

uhab3 = uhab_tmp |> 
  filter(str_detect(location, "Veg")) |> 
  filter(str_detect(trans, "D")) |> 
  ggplot() +
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = categ), alpha = 0.3) +
  geom_line(aes(x = datetime, .epred, linetype = categ), lineend = "round") +
  geom_point(aes(datetime, percent, color = categ)) +
  facet_rep_grid(cols = vars(species), rows = vars(habitat)) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_x_date("Month", date_minor_breaks = "month", labels = label_date_short()) +
  scale_y_continuous("Occurrence rate", limits = c(0, 100)) +
  ggtitle("(C) Vegetated: Deep") +
  theme_pubr() +
  scale_fill_viridis_d (end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  theme(
    plot.title       = element_text(size = rel(1.6)),
    legend.position  = "none",
    strip.background = element_rect(fill = "black"),
    strip.text.x     = element_text(face = "italic"),
    strip.text       = element_text(color = "white", size = rel(1.2)),
    axis.title       = element_blank(),
    axis.text        = element_blank()
  )

uhab4 = uhab_tmp |> 
  filter(str_detect(location, "Veg")) |> 
  filter(str_detect(trans, "S")) |> 
  ggplot() +
  geom_ribbon(aes(x = datetime, ymin = .lower, ymax = .upper, fill = categ), alpha = 0.3) +
  geom_line(aes(x = datetime, .epred, linetype = categ), lineend = "round") +
  geom_point(aes(datetime, percent, color = categ)) +
  facet_rep_grid(cols = vars(species), rows = vars(habitat)) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_x_date("Month", date_minor_breaks = "month", labels = label_date_short()) +
  scale_y_continuous("Occurrence rate", limits = c(0, 100)) +
  ggtitle("(D) Vegetated: Shallow") +
  theme_pubr() +
  scale_fill_viridis_d (end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  theme(
    plot.title       = element_text(size = rel(1.6)),
    legend.position  = "none",
    strip.background = element_rect(fill = "black"),
    strip.text.x     = element_blank(),
    strip.text       = element_text(color = "white", size = rel(1.2)),
    axis.text.y      = element_blank(),
    axis.title.y     = element_blank(),
    axis.text.x      = element_text(size = rel(1)),
    axis.title.x     = element_text(size = rel(1.5)),
    panel.spacing.y  = unit(-1.5, "lines")
  )

uhab1/uhab2|uhab3/uhab4


w = 400
h = 300

ggsave("Fig_5.pdf", width = w, height = h, units = "mm")

uhab_null = loo(uhab_mod_null_out)
uhab_mod  = loo(uhab_mod_out)

# Figure 6(Displacement and group GLM) -----------------------------------------

# Urchin displacement GLM (load .rds model fits from urchin_ecosystem_analysis.R)
disp_null_out = read_rds("tag_distance_null.rds")
disp_mod_out  = read_rds("tag_distance.rds")

obs_dist = df_tag |> 
  mutate(datetime = floor_date(datetime, "month")) |> 
  select(survey, location, trial, species, distance) |> 
  mutate(survey = factor(survey, levels = c("start", "mid", "end"), 
                         labels = c("0-hr", "12-hr", "24-hr"))) |> 
  mutate(species = recode(species, dset = "D. setosum", hcra = "H. crassispina"))

tag1 = 
  disp_mod_out$data |> add_epred_draws(disp_mod_out, seed = 2022) |> 
  group_by(species, location, survey) |> 
  mean_hdci(.epred) |> 
  mutate(species = recode(species, dset = "D. setosum", hcra = "H. crassispina")) |> 
  ggplot(aes(survey, .epred, color = species)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper), position = position_dodge(1)) +
  geom_point(aes(survey, distance),
             data = obs_dist, 
             alpha = 0.2,
             position = position_jitter(width = 0.1)) +
  facet_rep_wrap("location") +
  scale_y_continuous("Linear displacement (m)") +
  xlab("Survey") +
  ggtitle("(A)") +
  theme_pubr() +
  scale_color_viridis_d(end = 0.7,
                        labels = c(expression(italic("D. setosum")),
                                   expression(italic("H. crassispina")))) +
  theme(
    strip.background  = element_rect(fill = "black"),
    strip.text        = element_text(colour = "white", size = 13),
    legend.title      = element_blank(),
    legend.background = element_blank(),
    legend.position   = c(0.1, 0.9),
    legend.spacing    = unit(-2, "mm"),
    legend.key.size   = unit(0.45, "lines"),
    legend.text.align = 0,
    axis.text.x       = element_blank(),
    axis.title.x      = element_blank(),
    axis.text         = element_text(size = 10)
  )

# Urchin group composition GLM (load .rds model fits from urchin_ecosystem_analysis.R)
group_null_out  = read_rds("tag_group_null.rds")
group_mod_out2  = read_rds("tag_group.rds")

tag2 = 
  group_mod_out$data |> add_epred_draws(group_mod_out, seed = 2022) |> 
  group_by(location, survey) |> 
  mean_hdci(.epred) |> 
  ggplot(aes(survey, .epred)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper, color = "#330033"), 
                  position = position_nudge(x = -0.25)) +
  geom_point(aes(survey, individuals, color = "#330033"), data = df_tag |> drop_na() |> 
               filter(str_detect(species, "dset")) |> 
               mutate(survey = factor(survey, levels = c("start", "mid", "end"), 
                                      labels = c("0-hr", "12-hr", "24-hr"))),
             alpha = 0.2,
             position = position_jitter(width = 0.1)) +
  facet_rep_grid(cols = vars(location)) +
  scale_y_continuous("Group composition (indiv.)", limits = c(0, 25)) +
  scale_color_manual(values = "#330033") +
  xlab("Survey") +
  ggtitle("(B)") +
  theme_pubr() + 
  theme(
    legend.position  = "none",
    strip.text.x     = element_blank(),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(color = "white"),
    axis.text        = element_text(size = 10)
  )

tag1/tag2

w = 200
h = 200
ggsave("Fig_6.pdf", width = w, height = h, units = "mm")

disp_null  = loo(disp_null_out)
disp_mod   = loo(disp_mod_out)  
group_null = loo(group_null_out)
group_mod  = loo(group_mod_out)

# Figure 7(Urchin thresholds) --------------------------------------------------
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

u_sqm |> 
  dplyr::select(datetime, location, trans, total_per_sqm) |> 
  unique() |> 
  ggplot(aes(datetime, total_per_sqm)) +
  geom_point(size = 0.5) +
  geom_line() +
  facet_rep_grid(cols = vars(location), rows = vars(trans)) +
  theme_pubr() +
  geom_hline(yintercept = 8, linetype = "dashed", color = "firebrick",   size = 0.2) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "darkorange",  size = 0.2) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "chartreuse4", size = 0.2) +
  scale_y_continuous(expression(paste("Total urchins " ~m^{-2}))) +
  scale_x_datetime("Month", date_minor_breaks = "month", labels = label_date_short()) +
  theme(
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(color = "white"),
    axis.text        = element_text(size = 10)
  )

w = 150
h = 150

ggsave(path = "folder_figures", "Fig_7.pdf", width = w, height = h, units = "mm")

# Figure 8(Load old and new Arikawa satellite photos) --------------------------
# screen shots taken from: https://www.gsi.go.jp/johofukyu/johofukyu210322.html
# Screen shots are in the mapping_files folder together with other mapping resources
pic3 = "arikawa_old.png" |> image_read()
pic4 = "arikawa_new.png" |> image_read()

A3  = image_resize(pic3,  geometry = "500x500!")
A4  = image_resize(pic4,  geometry = "500x500!")

text1 = "(A)"
text2 = "°"
text3 = "(B)"

# Text background
w   = image_info(A3) |> pull(width)/8 # background width
h   = w # background height
w1  = image_info(A3) |> pull(width)
h1  = image_info(A3) |> pull(height)
bg2 = image_blank(w, h,   color = "rgba(0, 0, 0, 0)") # 100% transparency
bg3 = image_blank(w1, h1, color = "rgba(0, 0, 0, 0)")

bg_in_1 = bg2 |> image_annotate(text   = text1, gravity = "center", location = "+0+0",
                                size   = 30,       font = "Noto Serif",
                                weight = 500,     color = "black")
bg_in_1.1 = bg3 |> image_annotate(text   = text2, gravity = "northwest", location = "+0+0",
                                  size   = 30,       font = "Noto Serif",
                                  weight = 500,     color = "red")
image3 = image_composite(A3, bg_in_1, operator = "atop", offset = "+0+0") |> 
  image_border(color = "black", geometry = "1x1")
image3.1 = image_composite(image3, bg_in_1.1,   operator = "atop", offset = "+230+140")
image4.1 = image_composite(image3.1, bg_in_1.1, operator = "atop", offset = "+160+140")
image5.1 = image_composite(image4.1, bg_in_1.1, operator = "atop", offset = "+140+250")
image6.1 = image_composite(image5.1, bg_in_1.1, operator = "atop", offset = "+150+295")
image7.1 = image_composite(image6.1, bg_in_1.1, operator = "atop", offset = "+260+260")
image8.1 = image_composite(image7.1, bg_in_1.1, operator = "atop", offset = "+170+210")


bg_in_2 = bg2 |> image_annotate(text = text3, gravity = "center", location = "+0+0",
                                size = 30, font = "Noto Serif",
                                weight = 500, color = "black")
image4 = image_composite(A4, bg_in_2, operator = "atop", offset = "+0+0") |> 
  image_border(color = "black", geometry = "1x1")


arikawa_sat = image_append(c(image8.1, image4), stack = F)

image_write(arikawa_sat, format = "pdf", path = "./Fig_8.pdf", density = 100)

# Figure 9(Arikawa average winter temperature) ---------------------------------
temp_mean |> 
  filter(year != 2022) |> 
  mutate(month = as.numeric(month),
         col = ifelse(month<=3 & month>=1, "ok", "ng")) |> 
  filter(str_detect(col, "ok")) |> 
  mutate(datetime = floor_date(datetime, "month")) |>
  filter(!str_detect(position, "surface")) |> 
  group_by(year) |> 
  summarise(mean = mean(mean), sd = mean(sd)) |>
  ggplot(aes(year, mean, group = 1)) +
  geom_pointrange(aes(ymin = mean-sd, ymax = mean+sd)) +
  geom_line() +
  scale_y_continuous("Temperature (°C)") +
  xlab("Year") +
  theme_pubr() 

w = 150
h = 150

ggsave("Fig_9.pdf", w = w, h = h, units = "mm")  

################################# TABLES #######################################
# Table 1(Wave GAM) ------------------------------------------------------------
table_1_wave = wave_out |> 
  mutate(month = month(datetime, label = TRUE),
         year = year(datetime),
         Month = paste(year, month, sep = "-")) |> 
  dplyr::select(Month, Location = location, "Mean wave height" = .epred, Lower = .lower, Upper = .upper) |> 
  mutate_if(is.numeric, ~format(round(., 2))) |> 
  arrange(Location) 

write_csv(table_1_wave, file = "table_1_wave.csv")

# Table 2(Light GAM) -----------------------------------------------------------
table_2_light = light_out |> 
  mutate(month = month(datetime, label = TRUE),
         year = year(datetime),
         Month = paste(year, month, sep = "-")) |> 
  dplyr::select(Month, Location = location, Transect = trans, `Mean PPFD` = .epred, Lower = .lower, Upper = .upper) |> 
  mutate_if(is.numeric, ~format(round(., 2))) |> 
  arrange(Location, Transect)

write_csv(table_2_light, file = "table_2_light.csv")

# Table 3(Temperature GAM) -----------------------------------------------------
table_3_temp = temp_out |> 
  mutate(month = month(datetime, label = TRUE),
         year = year(datetime),
         Month = paste(year, month, sep = "-")) |> 
  dplyr::select(Month, Location = location, Transect = trans, `Mean temperature` = .epred, Lower = .lower, Upper = .upper) |> 
  mutate_if(is.numeric, ~format(round(., 2))) |> 
  arrange(Location, Transect)

write_csv(table_3_temp, file = "table_3_temp.csv")

# Table 4(ELPD LOO) ------------------------------------------------------------  
env1 = loo_compare(wave_mod, wave_null)   |> as_tibble() |> mutate(model = "Wave height GAM") |> slice(2)
env2 = loo_compare(light_mod, light_null) |> as_tibble() |> mutate(model = "Light (PPFD) GAM") |> slice(2)  
env3 = loo_compare(temp_mod, temp_null)   |> as_tibble() |> mutate(model = "Temperature GAM") |> slice(2) 

rug  = loo_compare(rug_mod, rug_null) |> as_tibble() |> mutate(model = "Rugosity GLM") |> slice(2)  

benth1 = loo_compare(benth_c_mod, benth_c_null) |> as_tibble() |> mutate(model = "Coralline GAM") |> slice(2)
benth2 = loo_compare(benth_m_mod, benth_m_null) |> as_tibble() |> mutate(model = "Macroalgal GAM") |> slice(2)
benth3 = loo_compare(benth_s_mod, benth_s_null) |> as_tibble() |> mutate(model = "Substrate GAM") |> slice(2)
benth4 = loo_compare(benth_t_mod, benth_t_null) |> as_tibble() |> mutate(model = "Turf GAM") |> slice(2)

sqm1 = loo_compare(sqm_dsav_mod, sqm_dsav_null) |> as_tibble() |> mutate(model = "D. savignyi density GAM") |> slice(2)
sqm2 = loo_compare(sqm_dset_mod, sqm_dset_null) |> as_tibble() |> mutate(model = "D. setosum density GAM") |> slice(2)
sqm3 = loo_compare(sqm_hcra_mod, sqm_hcra_null) |> as_tibble() |> mutate(model = "H. crassispina density GAM") |> slice(2)

bio1 = loo_compare(bio_dsav_mod, bio_dsav_null) |> as_tibble() |> mutate(model = "D. savignyi biomass GAM") |> slice(2)
bio2 = loo_compare(bio_dset_mod, bio_dset_null) |> as_tibble() |> mutate(model = "D. setosum biomass GAM") |> slice(2)
bio3 = loo_compare(bio_hcra_mod, bio_hcra_null) |> as_tibble() |> mutate(model = "H. crassispina biomass GAM") |> slice(2)

hab   = loo_compare(uhab_mod, uhab_null)   |> as_tibble() |> mutate(model = "Microhabitat GAM") |> slice(2)
disp  = loo_compare(disp_mod, disp_null)   |> as_tibble() |> mutate(model = "Linear displacement GLM") |> slice(2)
group = loo_compare(group_mod, group_null) |> as_tibble() |> mutate(model = "Group-size GLM") |> slice(2)


table_4 = bind_rows(env1, env2, env3, rug, benth1, benth2, benth3, benth4, 
                     sqm1, sqm2, sqm3, bio1, bio2, bio3, hab, disp, group) |> 
  dplyr::select(Model = model, "ELPD difference" = elpd_diff, "Standard error" = se_diff) |> 
  mutate_at(c("ELPD difference", "Standard error"), ~format(round(., 2)))

write_csv(table_4, "table_4_elpdloo.csv")


# Table 5(Benthic cover GAM) ---------------------------------------------------
table_5_cover = benth_out |> 
  mutate(month = month(datetime, label = TRUE),
         year = year(datetime),
         Month = paste(year, month, sep = "-")) |> 
  dplyr::select(Month, Location = location, Transect = trans, 
                Element = element, `Percent cover` = .epred, Lower = .lower, Upper = .upper) |> 
  mutate_if(is.numeric, ~format(round(., 2))) |> 
  arrange(Location, Transect)

write_csv(table_5_cover, file = "table_5_cover.csv")

# Table 6(Urchin density GAM) --------------------------------------------------
table_6_density = urchin_count |> 
  mutate(month = month(datetime, label = TRUE),
         year = year(datetime),
         Month = paste(year, month, sep = "-")) |> 
  dplyr::select(Month, Location = location, Transect = trans, Species = species,
                Density = .epred, Lower = .lower, Upper = .upper) |> 
  mutate_if(is.numeric, ~format(round(., 2))) |> 
  arrange(Location, Transect)

write_csv(table_6_density, file = "table_6_density.csv")

# Table 7(Urchin biomass GAM) -------------------------------------------------- 
table_7_biomass = urchin_biomass |> 
  mutate(month = month(datetime, label = TRUE),
         year = year(datetime),
         Month = paste(year, month, sep = "-")) |> 
  dplyr::select(Month, Location = location, Transect = trans, Species = species, Biomass = .epred,
                Lower = .lower, Upper = .upper) |> 
  mutate_if(is.numeric, ~format(round(., 2))) |> 
  arrange(Location, Transect)

write_csv(table_7_biomass, file = "table_7_biomass.csv")

# Table 8(Microhabitat GAM) ----------------------------------------------------
table_8_microhabitat = uhab_tmp |> 
  mutate(month = month(datetime, label = TRUE),
         year = year(datetime),
         Month = paste(year, month, sep = "-")) |> 
  dplyr::select(datetime, Month, Location = location, Transect = trans, 
                Species = species, `Size class` = categ, Microhabitat = habitat, 
                `Occurrence rate` = .epred, Lower = .lower, Upper = .upper) |> 
  mutate_if(is.numeric, ~format(round(., 2))) |> 
  arrange(datetime, Location, Transect) |> 
  dplyr::select(-datetime)

write_csv(table_8_microhabitat, file = "table_8_microhabitat.csv")

# Table 9(Linear displacement) -------------------------------------------------
table_9_tagging = 
  disp_mod_out$data |> add_epred_draws(disp_mod_out, seed = 2022) |>
  group_by(species, location, survey) |>
  mean_hdci(.epred) |>
  dplyr::select(Location = location, Species = species, Survey = survey, 
                Displacement = .epred, Lower = .lower, Upper = .upper) |> 
  mutate(Species = recode(Species, dset = "D. setosum", hcra = "H. crassispina")) |> 
  arrange(Location) |> 
  mutate_if(is.numeric, ~ format(round(., 2)))

write_csv(table_9_tagging, "table_9_tag_linear.csv")

# Table 10(group size) ---------------------------------------------------------
table_10_tagging = 
  group_mod_out$data |> add_epred_draws(group_mod_out, seed = 2022) |> 
  group_by(location, survey) |> 
  mean_hdci(.epred) |> 
  mutate(species = rep_len("D. setosum", length.out = 6)) |> 
  dplyr::select(Location = location, Species = species, Survey = survey, 
                Individuals = .epred, Lower = .lower, Upper = .upper) |> 
  arrange(Location) |> 
  mutate_if(is.numeric, ~ format(round(., 2)))

write_csv(table_10_tagging, "table_10_tag_group.csv")

# Table 11(Recapture rate) -----------------------------------------------------
table_11_recapture = df_tag |> 
  filter(str_detect(remarks, "ok")) |> 
  group_by(species, location, survey) |>
  mutate(survey = factor(survey, 
                         levels = c("start", "mid", "end"), 
                         labels = c("0-hr", "12-hr", "24-hr")),
         species = recode(species, dset = "D. setosum", hcra = "H. crassispina")) |> 
  distinct(trial) |> 
  summarise(n = length(species), .groups = "drop") |> 
  mutate(`Locate rate (%)` = (n/20)*100,
         `Loss rate (%)` = 100 - `Locate rate (%)`) |> 
  arrange(location) |> 
  dplyr::select(Location = location, Species = species, Survey = survey, n,
                `Locate rate (%)`, `Loss rate (%)`)
  
write_csv(table_11_recapture, file = "table_11_recapture.csv")

######################### SUPPLEMENTARY FIGURES ################################
# Figure S1(Wind fetch) --------------------------------------------------------
# For details see: https://www.frontiersin.org/articles/10.3389/fmars.2022.861932/full#B88 
calc_intersection = function(origin, map_layer, fetch_limits) {
  # This function creates a line from the origin to the fetch limit.
  # Then it will determine the points where the line and the map polygons
  # intersects. There can be more than one intersection if the line 
  # has sufficient length. The trick is to determine the point closest to the
  # origin so that it can be used as the end point for the fetch. If
  # there is no intersection, then the end point for the fetch is the 
  # fetch limit.
  fl = fetch_limits[1:2]
  or = st_coordinates(origin)
  X = matrix(c(or, fl), ncol = 2, byrow = T)
  lst = st_linestring(X) |> st_sfc(crs = st_crs(map_layer))
  int = st_intersection(lst, map_layer) 
  
  xy  = st_coordinates(int) |> as.matrix()
  distance = function(x1, x2) {
    sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2 )
  }
  # If there are not intersections, then the xy matrix has dimensions [0,2]
  j = dim(xy)
  if(j[1] > 0) {
    dis = apply(xy, 1, distance, x2 = fl)
    max_dis = max(dis)
    n = which(near(dis, max_dis)) 
    xy = xy[n,]
    X = matrix(c(or, xy[1:2]), ncol = 2, byrow =T)
    lst = st_linestring(X) |> st_sfc(crs = st_crs(map_layer))
    len = st_length(lst) 
    z = lst |> st_as_sf() |> mutate(length = len)
  } else {
    X = matrix(c(or, fl), ncol = 2, byrow =T)
    lst = st_linestring(X) |> st_sfc(crs = st_crs(map_layer))
    len = st_length(lst) 
    z = lst |> st_as_sf() |> mutate(length = len) 
  }
  z
}

calc_circle = function(map_layer, max_dist=30, n_vectors = 9) {
  # Calculate the fetch limits.
  delta_theta =  360 / (n_vectors * 4)
  theta = seq(0, 360, by = delta_theta) |> head(n = -1)
  max_dist = max_dist * 1000 # Convert from km to m.
  max_dist = units::set_units(max_dist, "m")
  d_bff = st_buffer(map_layer, dist = max_dist, nQuadSegs = n_vectors) 
  fetch_ends = st_coordinates(d_bff) |> head(n = -1) 
  list(d_bff = d_bff, fetch_limits = fetch_ends[order(theta), ], theta = theta)
}
# GPS coordinates to determine wind fetch 
vegetatedgps = matrix(rev(c(32.987887, 129.118768)), ncol = 2)
isoyakegps    = matrix(rev(c(32.9859638, 129.119355)), ncol = 2)

gps_info = 
  rbind(vegetatedgps, isoyakegps) |> 
  as_tibble(.name_repair = ~c("long", "lat")) |> 
  mutate(name = c("vegetatedgps", "isoyakegps")) |> 
  mutate(label = str_to_sentence(str_remove(name, pattern = "(gps)"))) 

gps_info = gps_info |> 
  mutate(label2 = str_to_sentence(label)) |> 
  mutate(label2 = str_glue("{label2} {ifelse(str_detect(label2, 'Bise'), 'Point', 'Bay')}"))

# Prepare the Coordinate Reference System to be EPSG:4326 (Which is WGS 84)
# See st_crs(4326) for details
gps_info = gps_info |> dplyr::select(long, lat, name) |> st_as_sf(coords = c("long", "lat"), crs = 4326, agr = "constant")

# Load the map shape files 
# The map uses the ITRF94 system (st_crs(map_poly))
map_poly = read_sf("~/Lab_Data/Japan_map_data/GSI/polbnda_jpn.shp")
map_poly = map_poly |> dplyr::select(nam, geometry)

# Convert the CRS to EPSG:2450 
map_poly = st_transform(map_poly, st_crs(2450))
gps_info  = st_transform(gps_info, st_crs(2450))

# Do the analysis one location at a time. 
ptsize = 1
max_dist = 10 # In km
n_vectors = 3*9 # The number of vectors in every quadrant.


location = "Vegetated habitat"
polygon_layer = subset(map_poly, str_detect(nam, "Nagasaki")) |> st_union() 
site_layer    = subset(gps_info, str_detect(name, "veg"))
fetch_limits = calc_circle(site_layer, max_dist = max_dist, n_vectors = n_vectors)
fout = fetch_limits$fetch_limits |> as_tibble() |> 
  mutate(fe  = map2(X,Y,function(x,y) cbind(x,y))) |> 
  mutate(geometry = purrr::map(fe, calc_intersection, origin = site_layer, map_layer = polygon_layer))
fout = fout |> select(geometry) |>  unnest(geometry) |> st_as_sf()
temp_layer = st_crop(polygon_layer, st_bbox(fetch_limits$d_bff))
mean_fetch = fout |> pull(length) |> mean() |> as.numeric()
sd_fetch = fout |> pull(length) |> sd() |> as.numeric()

max_fetch = fout |> pull(length) |> as.numeric()
man_n = sum(near(max_fetch, max_dist * 1000))
tot_n = length(max_fetch)

p1 = ggplot() + 
  geom_sf(data = temp_layer, color = NA) +
  geom_sf(data = fout) +
  geom_sf(data = site_layer, color = "red", size = ptsize) +
  labs(title = str_glue("The mean ± sd fetch for {location} is {format(mean_fetch, digits = 4)} ± {format(sd_fetch, digits = 4)} m."),
       subtitle = str_glue("{man_n} out of {tot_n} vectors were at the upper limit.")) +
  theme_pubr()

location = "Isoyake habitat"
polygon_layer = subset(map_poly, str_detect(nam, "Nagasaki")) |> st_union() 
site_layer    = subset(gps_info, str_detect(name, "iso"))
fetch_limits = calc_circle(site_layer, max_dist = max_dist, n_vectors = n_vectors)
fout = fetch_limits$fetch_limits |> as_tibble() |> 
  mutate(fe  = map2(X,Y,function(x,y) cbind(x,y))) |> 
  mutate(geometry = map(fe, calc_intersection, origin = site_layer, map_layer = polygon_layer))
fout = fout |> select(geometry) |>  unnest(geometry) |> st_as_sf()
temp_layer = st_crop(polygon_layer, st_bbox(fetch_limits$d_bff))
mean_fetch = fout |> pull(length) |> mean() |> as.numeric()
sd_fetch = fout |> pull(length) |> sd() |> as.numeric()

max_fetch = fout |> pull(length) |> as.numeric()
man_n = sum(near(max_fetch, max_dist * 1000))
tot_n = length(max_fetch)

p2 = ggplot() + 
  geom_sf(data = temp_layer, color = NA) +
  geom_sf(data = fout) +
  geom_sf(data = site_layer, color = "red", size = ptsize) +
  labs(title = str_glue("The mean ± sd fetch for {location} is {format(mean_fetch, digits = 4)} ± {format(sd_fetch, digits = 4)} m."),
       subtitle = str_glue("{man_n} out of {tot_n} vectors were at the upper limit.")) +
  theme_pubr()

p1 + p2 + plot_layout(ncol = 2)

ggsave("Fig_S1.pdf", width = 4*80, height = 3*80, units = "mm")

# Figure S2(Environmental GAM kernel density) ----------------------------------
pp_wave  = pp_check(m_wave_out, type = "hist",ndraws = 5) +
  ggtitle("(A) Wave height") +
  theme(
    legend.position   = c(0.2, 0.93),
    legend.key.size   = unit(0.5, "lines"),
    legend.direction  = "horizontal",
    legend.background = element_blank(),
    panel.border      = element_rect(colour = "black", fill=NA),
    legend.text       = element_text(size = 8), 
  )
pp_light = pp_check(m_light_out, type = "hist",ndraws = 5) +
  ggtitle("(B) Light (PPFD)") +
  theme(
    legend.position = "none",
    panel.border    = element_rect(colour = "black", fill=NA),
    legend.text     = element_text(size = 15), 
  )
pp_temp  = pp_check(m_temp_out, type = "hist", ndraws = 5) +
  ggtitle("(C) Temperature") + 
  theme(
    panel.border    = element_rect(colour = "black", fill=NA),
    legend.position = "none", 
    text            = element_text(size = 10),
  )

pp_wave/pp_light/pp_temp
w = 200
h = 180

ggsave("Fig_S2.pdf", width = w, height = h, units = "mm")

# Figure S3(Environmental GAM MCMC rank plots) ---------------------------------
mcmc_wave =
  mcmc_rank_hist(m_wave_out, pars = vars("b_Intercept", "b_locationVegetated")) +
  ggtitle("(A) Wave") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )

mcmc_light =
  mcmc_rank_hist(m_light_out, pars = vars("b_Intercept", "b_locationVegetated", "b_transShallow")) + 
  ggtitle("(B) Light") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )

mcmc_temp = 
  mcmc_rank_hist(m_temp_out, pars = vars("b_Intercept", "b_locationVegetated", "b_transShallow")) + 
  ggtitle("(C) Temperature") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 8),
  )

mcmc_wave/mcmc_light/mcmc_temp

w = 250
h = 300

ggsave("Fig_S3.pdf",width = w, height = h, units = "mm")

# Figure S4(Benthic rugosity GLM kernel density) -------------------------------
pp_check(m_rug_out, type = "hist", ndraws = 5) +
  ggtitle("Rugosity GLM") +
  theme(
    legend.position   = c(0.25, 0.95),
    legend.direction  = "horizontal",
    legend.key.size   = unit(0.5, "lines"),
    legend.text       = element_text(size   = 8),
    legend.background = element_blank(),
    panel.border      = element_rect(colour = "black", fill=NA),
    text              = element_text(size   = 8)
  )
w = 200
h = 100

ggsave("Fig_S4.pdf", width = w, height = h, units = "mm")

# Figure S5(Benthic rugosity GLM MCMC rank plot) -------------------------------
mcmc_rank_hist(m_rug_out, pars = vars("b_Intercept", "b_locationVegetated", "b_transshallow")) +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 7),
  )

w = 200
h = 130

ggsave("Fig_S5.pdf", width = w, height = h, units = "mm")

# Figure S6(Benthic cover GAM kernel density plots) ----------------------------
pp_cor = pp_check(m_benth_c_out, type = "hist", ndraws = 5) + 
  ggtitle("(A) Coralline") +
  theme(
    legend.position   = c(0.13, 0.92),
    legend.direction  = "horizontal",
    legend.background = element_blank(),
    legend.key.size   = unit(0.5, "lines"),
    panel.border      = element_rect(colour = "black", fill=NA),
    legend.text       = element_text(size = 8), 
    text              = element_text(size = 10)
  )
pp_mac = pp_check(m_benth_m_out, type = "hist", ndraws = 5) + 
  ggtitle("(B) Macroalgae") + 
  theme(
    legend.position = "none", 
    panel.border    = element_rect(colour = "black", fill=NA),
    text            = element_text(size = 10),
  )
pp_sub = pp_check(m_benth_s_out, type = "hist", ndraws = 5) + 
  ggtitle("(C) Substrate") + 
  theme(
    panel.border    = element_rect(colour = "black", fill=NA),
    legend.position =   "none", 
    text            = element_text(size = 10)
  )
pp_tur = pp_check(m_benth_t_out, type = "hist", ndraws = 5) +
  ggtitle("(D) Turf algae") + 
  theme(
    panel.border    = element_rect(colour = "black", fill=NA),
    legend.position = "none", 
    text            = element_text(size = 10)
  )

pp_cor/pp_mac|pp_sub/pp_tur

w = 250
h = 180

ggsave("Fig_S6.pdf", width = w, height = h, units = "mm")

# Figure S7(Benthic cover GAM MCMC rank plots) ---------------------------------
rank_cor = mcmc_rank_hist(m_benth_c_out, pars = vars("b_Intercept", "b_locationVegetated")) +
  ggtitle("(A) Coralline algae") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )
rank_mac = mcmc_rank_hist(m_benth_m_out, pars = vars("b_Intercept", "b_locationVegetated")) +
  ggtitle("(B) Macroalgae") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )
rank_sub = mcmc_rank_hist(m_benth_s_out, pars = vars("b_Intercept", "b_locationVegetated")) +
  ggtitle("(C) Substrate") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )
rank_tur = mcmc_rank_hist(m_benth_t_out, pars = vars("b_Intercept", "b_locationVegetated")) +
  ggtitle("(D) Turf algae") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
  )

rank_cor/rank_mac/rank_sub/rank_tur

w = 180
h = 280

ggsave("Fig_S7.pdf", width = w, height = h, units = "mm")

# Figure S8(Urchin density and biomass GAM kernel density plots) ---------------
pp_sqm_dsav = pp_check(dsav_sqm_out, type = "hist", ndraws = 5) +
  ggtitle("(A) D. savignyi density") + 
  theme(
    legend.position   = c(0.2, 0.92),
    legend.key.size   = unit(0.5, "lines"),
    legend.direction  = "horizontal",
    legend.text       = element_text(size = 8),
    legend.background = element_blank(),
    panel.border      = element_rect(colour = "black", fill=NA),
    text              = element_text(size = 10)
  )
pp_sqm_dset = pp_check(dset_sqm_out, type = "hist", ndraws = 5) +
  ggtitle("(B) D. setosum density") +
  theme(
    panel.border    = element_rect(colour = "black", fill=NA),
    legend.position = "none", 
    text            = element_text(size = 10)
  )
pp_sqm_hcra = pp_check(hcra_sqm_out, type = "hist", ndraws = 5) +
  ggtitle("(C) H. crassispina density") + 
  theme(
    panel.border    = element_rect(colour = "black", fill=NA),
    legend.position = "none", 
    text            = element_text(size = 10)
  )

pp_bio_dsav = pp_check(dsav_bio_out, type = "hist", ndraws = 5) +
  ggtitle("(D) D. savignyi biomass") + 
  theme(
    legend.position   = "none",
    legend.background = element_blank(),
    panel.border      = element_rect(colour = "black", fill=NA),
    text              = element_text(size = 10)
  )
pp_bio_dset = pp_check(dset_bio_out, type = "hist", ndraws = 5) +
  ggtitle("(E) D. setosum biomass") +
  theme(
    panel.border    = element_rect(colour = "black", fill=NA),
    legend.position = "none", 
    text            = element_text(size = 10)
  )
pp_bio_hcra = pp_check(hcra_bio_out, type = "hist", ndraws = 5) +
  ggtitle("(F) H. crassispina biomass") +
  theme(
    panel.border    = element_rect(colour = "black", fill=NA),
    legend.position = "none", 
    text            = element_text(size = 10)
  )

w = 250
h = 200

pp_sqm_dsav/pp_sqm_dset/pp_sqm_hcra | pp_bio_dsav/pp_bio_dset/pp_bio_hcra

ggsave("Fig_S8.pdf", width = w, height = h, units = "mm")

# Figure S9(Urchin density GAM MCMC rank plots) --------------------------------
mcmc_sqm_dsav = 
  mcmc_rank_hist(dsav_sqm_out, pars = vars("b_Intercept", "b_locationVegetated", "b_transShallow")) +
  ggtitle("(A) D. savignyi density") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
  )

mcmc_sqm_dset = 
  mcmc_rank_hist(dset_sqm_out, pars = vars("b_Intercept", "b_locationVegetated", "b_transShallow")) +
  ggtitle("(B) D. setosum density") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
  )

mcmc_sqm_hcra =
  mcmc_rank_hist(hcra_sqm_out, pars = vars("b_Intercept", "b_locationVegetated", "b_transShallow")) + 
  ggtitle("(C) H. crassispina density") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
  )

mcmc_sqm_dsav/mcmc_sqm_dset/mcmc_sqm_hcra

w = 180
h = 280

ggsave("Fig_S9.pdf", width = w, height = h, units = "mm")

# Figure S10(Urchin biomass GAM MCMC rank plots) --------------------------------
mcmc_bio_dsav = 
  mcmc_rank_hist(dsav_bio_out, pars = vars("b_Intercept", "b_locationVegetated", "b_transShallow")) +
  ggtitle("(A) D. savignyi biomass") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
  )

mcmc_bio_dset = 
  mcmc_rank_hist(dset_bio_out, pars = vars("b_Intercept", "b_locationVegetated", "b_transShallow")) +
  ggtitle("(B) D. setosum biomass") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
  )

mcmc_bio_hcra =
  mcmc_rank_hist(hcra_bio_out, pars = vars("b_Intercept", "b_locationVegetated", "b_transShallow")) + 
  ggtitle("(C) H. crassispina biomass") +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 6),
  )

mcmc_bio_dsav/mcmc_bio_dset/mcmc_bio_hcra

ggsave("Fig_S10.pdf", width = w, height = h, units = "mm")

# Figure S11(Microhabitat GAM kernel density plots) ----------------------------
pp_check(uhab_mod_out, type = "hist", ndraws = 5) +
  theme(
    legend.position   = c(0.13, 0.94),
    legend.key.size   = unit(0.5, "lines"),
    legend.direction  = "horizontal",
    legend.text       = element_text(size   = 8),
    legend.background = element_blank(),
    panel.border      = element_rect(colour = "black", fill=NA),
    text              = element_text(size   = 10)
  )

w = 200
h = 100

ggsave("Fig_S11.pdf", width = w, height = h, units = "mm")

# Figure S12(Microhabitat GAM MCMC rank plot) ----------------------------------

mcmc_rank_hist(uhab_mod_out, pars = vars("b_Intercept", "b_speciesdset", "b_specieshcra")) +
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 7)
  )
w = 200
h = 100

ggsave("Fig_S12.pdf", width = w, height = h, units = "mm")

# Figure S13(Linear displacement GLM kernel density plot) ----------------------
pp_check(disp_mod_out, type = "hist", ndraws = 5) +
  theme(
    panel.border      = element_rect(colour = "black", fill=NA),
    legend.direction  = "horizontal",
    legend.background = element_blank(),
    legend.key.size   = unit(0.5, "lines"), 
    legend.text       = element_text(size   = 8),
    legend.position   = c(0.13, 0.94),
    text              = element_text(size = 10)
  )

w = 200
h = 130

ggsave("Fig_S13.pdf", width = w, height = h, units = "mm")

# Figure S14(Linear displacement GLM MCMC rank plot) ---------------------------
mcmc_rank_hist(disp_mod_out, pars = vars("b_Intercept", "b_specieshcra")) + 
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
  )

w = 200
h = 100

ggsave("Fig_S14.pdf", width = w, height = h, units = "mm")

# Figure S15(Group-size GLM kernel density plot) -------------------------------
pp_check(group_mod_out, type = "hist", ndraws = 5) +
  theme(
    panel.border    = element_rect(colour = "black", fill=NA),
    legend.direction  = "horizontal",
    legend.background = element_blank(),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.5, "lines"),
    legend.position   = c(0.14, 0.94),
    text              = element_text(size = 10)
  )
w = 200
h = 130

ggsave("Fig_S15.pdf", width = w, height = h, units = "mm")

# Figure S16(Group-size GLM MCMC rank plot) ------------------------------------
mcmc_rank_hist(group_mod_out, pars = vars("b_Intercept", "b_locationVegetated")) + 
  theme(
    text         = element_text(size = 10),
    strip.text.y = element_text(size = 7),
  )
w = 200
h = 100

ggsave("Fig_S16.pdf", width = w, height = h, units = "mm")

