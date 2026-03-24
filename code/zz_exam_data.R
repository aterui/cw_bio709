# exam for biostats II BIO 709 --------------------------------------------

### Reproducibility
set.seed(123)

### Load required packages
library(tidyverse)
library(glmmTMB)
library(piecewiseSEM)
library(mgcv)

## glmm

### Design: 10 lakes × 10 observations each
lake <- rep(letters[1:10], each = 10)
n <- length(lake)

### Create numeric index for lake (avoids repeated factor conversion)
lake_id <- as.numeric(factor(lake))
n_lake  <- length(unique(lake))

### Lake-level latent variables (random-effect-like drivers)
g1 <- rnorm(n_lake, mean = 10, sd = 2.5)   # influences hb and s
g2 <- runif(n_lake, min = 1, max = 10)     # influences prod

### Observation-level predictors
cond      <- runif(n, min = 1, max = 20)
substrate <- runif(n, min = 0, max = 7)

### Data-generating process (hierarchical structure)
## prod depends on lake-level g2 and condition
prod <- rnorm(
  n,
  mean = g2[lake_id] + 0.5 * cond,
  sd   = 1
)

## hb depends on lake-level g1 and productivity
mu <- g1[lake_id] + 0.2 * prod
hb <- rnorm(n, mean = mu)

## Poisson response depends on g1 and productivity
lambda <- exp(g1[lake_id] * 0.05 + 0.1 * prod)
s <- rpois(n, lambda = lambda)

### Assemble dataset
df_s <- tibble(
  s = s,
  hb = hb,
  prod = prod,
  substrate = substrate,
  cond = cond,
  lake = lake
)

### Visualization: relationship between prod and count response by lake
ggplot(df_s) +
  geom_point(aes(x = prod, y = s, color = lake)) +
  facet_wrap(~lake)

### Fit GLMM (Gaussian with log-transformed response)
m_glmm <- glmmTMB(
  log(hb) ~ prod + substrate + (1 | lake),
  data   = df_s,
  family = gaussian
)

summary(m_glmm)

### Save simulated dataset
saveRDS(df_s, "data_fmt/data_lake_invert.rds")


## psem

### Component models reflecting assumed causal structure
## hb is driven by prod and substrate
m1 <- glmmTMB(
  hb ~ prod + substrate + (1 | lake),
  data = df_s
)

## prod is driven by environmental condition
m2 <- glmmTMB(
  prod ~ cond + (1 | lake),
  data = df_s
)

## count response depends on productivity
m3 <- glmmTMB(
  s ~ prod + (1 | lake),
  data   = df_s,
  family = poisson
)

### Combine into piecewise structural equation model
fit_psem <- psem(m1, m2, m3)


## gam

### Reproducibility for GAM simulation
set.seed(123)

### Time index (daily for one year)
t <- 1:365

### Exponential decay modifier (e.g., seasonal decline)
p <- 0.995^(0:364)

### Simulate seasonal emergence patterns for two sites
## Site 1: longer seasonal cycle + lower noise
y1 <- (2 +
         sin(seq(-0.5 * pi, 7/2 * pi, length = 365)) +
         rnorm(365, sd = 0.3)) * p

## Site 2: slightly shifted cycle + higher noise
y2 <- (2.2 +
         sin(seq(-0.5 * pi, 6/2 * pi, length = 365)) +
         rnorm(365, sd = 0.4)) * p

### Assemble dataset in long format (required for GAM with factor smooths)
df_emg <- tibble(
  t  = t,
  s1 = y1,
  s2 = y2
) %>%
  pivot_longer(
    cols      = c(s1, s2),
    names_to  = "site",
    values_to = "emergence"
  ) %>%
  mutate(site = factor(site))

### Visualization: temporal patterns by site
ggplot(df_emg, aes(x = t, y = emergence, color = site)) +
  geom_point()

### Fit GAM with site-specific smooths
## s(t, by = site): separate smooth for each site
## + site: parametric difference in intercept
m_gam <- gam(
  emergence ~ s(t, by = site) + site,
  data = df_emg
)

summary(m_gam)

### Save dataset
saveRDS(df_emg, "data_fmt/data_insect_emergence.rds")

## Ordination
nutrient <- rnorm(nrow(trees), -50 + 0.8 * trees$Height, sd = 1)
saveRDS(nutrient, "data_fmt/nutrient.rds")

## time series

df_nile <- tibble(
  year = time(Nile),
  discharge = as.numeric(Nile)
)

df_sunspot <- tibble(
  year = time(sunspot.year),
  sunspots = as.numeric(sunspot.year)
)
