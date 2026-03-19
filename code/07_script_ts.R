#' DESCRIPTION:
#' Script for time-series

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               forecast,
               lterdatasampler,
               daymetr,
               glarma)

url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_ts_anormaly.csv"
(df_ts <- read_csv(url))

## draw a figure
g_base <- df_ts %>% 
  ggplot(aes(x = year,
             y = anormaly)) +
  geom_line() +
  geom_point() +
  labs(
    y = "Anormaly",
    x = "Year"
  )

## regression
m_ts <- lm(anormaly ~ year,
           data = df_ts)

summary(m_ts)

## figure
g_base +
  geom_abline(intercept = coef(m_ts)[1],
              slope = coef(m_ts)[2]) +
  theme_bw()

## random-walk
y <- NULL  
y[1] <- 0

for (i in 1:99) {
  y[i + 1] <- y[i] + rnorm(1, mean = 0, sd = 1)
}

tibble(y = y,
       x = 1:length(y)) %>% 
  ggplot(aes(x = x,
             y = y)) +
  geom_point() +
  geom_line()

## Lake Huron data
df_huron <- tibble(
  year = time(LakeHuron),                # Extracts the time component (years) from the LakeHuron ts object
  water_level = as.numeric(LakeHuron)    # Converts LakeHuron values to numeric (from ts class)
) %>% 
  arrange(year)  

# Plot Lake Huron time series with a linear trend
df_huron %>% 
  ggplot(aes(x = year, y = water_level)) +
  geom_point(alpha = 0.25) +       # Semi-transparent points
  geom_line(linetype = "dotted") + # Dotted line connecting points
  geom_smooth(method = "lm",       # Linear trend line
              color = "black",
              linewidth = 0.5) +
  theme_bw() +
  labs(x = "Year", y = "Water Level")

## Autoregressive model
m_ar1 <- Arima(
  df_huron$water_level,
  order = c(1, 0, 0)
)

## fitted values
df_huron_ar1 <- df_huron %>% 
  mutate(fit = fitted(m_ar1) %>% 
           as.numeric())

df_huron_ar1 %>% 
  ggplot() +
  geom_point(aes(x = year,
                 y = water_level),
             alpha = 0.25) +
  geom_line(aes(x = year,
                y = fit),
            color = "steelblue") +
  theme_bw()

## Moving-average model
m_ma1 <- Arima(
  df_huron$water_level,
  order = c(0, 0, 1)
)

## ARMA model
m_arma1 <- Arima(
  df_huron$water_level,
  order = c(1, 0, 1)
)

## ARIMA model
m_arima1 <- Arima(
  df_huron$water_level,
  order = c(1, 1, 0),
)

## model selection
auto.arima(
  y = df_huron$water_level,
  stepwise = FALSE,
  ic = "aic"
)

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Bison Body Mass, Climate, and Time-Series Analysis
# ============================================================

library(lterdatasampler)

# The "knz_bison" dataset contains long-term monitoring data
# on bison captured at Konza Prairie Biological Station.
#
# ------------------------------------------------------------
# Key columns may include:
# rec_year      : Year of capture
# animal_sex    : Sex of the individual (e.g., female, male)
# animal_weight : Body mass of bison
# ------------------------------------------------------------
#
# In this exercise, you will explore long-term trends in bison
# body mass and evaluate how climate variability may influence
# weight dynamics over time.

# 1. Explore the structure of the knz_bison dataset.
#    - Inspect variable types and missing values.
#    - Reformat variables as needed for analysis.

# 2. Subset the data to include observations from 1994–2012.

# 3. Calculate the average body mass for female and male bison
#    for each year in the selected time period.

# 4. Obtain climate data from the daymetr dataset.
#    - Identify relevant climate variables (e.g., temperature,
#      precipitation).
#    - Associate climate data with knz_bison by year.
#    - Coordinates: Lat 39.09300	Lon -96.57500

# 5. Perform a time-series analysis to examine whether selected
#    climate variables influence annual bison body mass.
#    - Consider temporal autocorrelation and lag effects.
#    - Model males and females separately

# 6. Using your fitted model, compare observed bison body mass
#    with predicted values for the period 2014–2020.
#    - Evaluate model performance and discuss sources of uncertainty.
