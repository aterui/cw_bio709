#' DESCRIPTION:
#' Script for Constrained Ordination

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               GGally,
               vegan)

data("varespec", "varechem")

# rename
m_y <- varespec
colnames(m_y) <- str_to_lower(colnames(m_y))

df_env <- as_tibble(varechem) %>% 
  janitor::clean_names()

# visualization
m_y %>% 
  ggpairs(
    progress = FALSE,
    column = 1:3,
    aes(alpha = 0.25)
  ) +
  theme_bw()

# perform RDA
(obj_rda <- rda(m_y ~ n + p + ca,
                data = df_env))

# statistical test
anova.cca(obj_rda,
          by = "margin",
          permutations = 999)

# visualization
# community data
df_rda <- scores(obj_rda,
                 display = "site",
                 scaling = 2) %>% 
  bind_cols(df_env) %>% 
  janitor::clean_names()

# env vectors
df_bp <- scores(obj_rda,
                display = "bp",
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

# Create a ggplot2 ordination plot
# - Points represent sites positioned by their constrained community composition
# - Color gradient reflects the nitrogen (n) concentration at each site
df_rda %>% 
  ggplot(aes(x = rda1,
             y = rda2)) +        # color sites by nitrogen level
  geom_point(aes(color = n)) +
  geom_segment(data = df_bp,
               aes(x = 0, xend = rda1 * 10, # 10 is arbitrary scaling for visualization
                   y = 0, yend = rda2 * 10),
               arrow = arrow(length = unit(0.2, "cm"))
  ) +
  geom_text(data = df_bp,
            aes(x = rda1 * 10.5,    # slightly beyond arrow tip
                y = rda2 * 10.5,
                label = variable),  # or use a variable column
            size = 4) +
  theme_bw() +
  labs(x = "RDA1",
       y = "RDA2",
       color = "Nitrogen") +
  scale_color_viridis_c()

# perform dbRDA
(obj_db <- dbrda(m_y ~ n + p + ca,
                 data = df_env,
                 distance = "bray"))

# test with anova.cca
anova.cca(obj_db,
          by = "margin",
          permutation = 999)

# visualization
df_db <- scores(obj_db,
                display = "sites",
                scaling = 2) %>% 
  as_tibble() %>% 
  bind_cols(df_env) %>% 
  janitor::clean_names()

df_bp <- scores(obj_db,
                display = "bp",
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

df_db %>% 
  ggplot(aes(x = db_rda1,
             y = db_rda2)) +
  geom_point(aes(color = n)) +
  geom_segment(data = df_bp,
               aes(x = 0, xend = db_rda1,
                   y = 0, yend = db_rda2),
               arrow = arrow(length = unit(0.2, "cm"))
  ) +
  geom_text(data = df_bp,
            aes(x = db_rda1 * 1.1,    # slightly beyond arrow tip
                y = db_rda2 * 1.1,
                label = variable),  # or use a variable column
            size = 4) +
  theme_bw() +
  labs(x = "dbRDA1",
       y = "dbRDA2",
       color = "Nitrogen") +
  scale_color_viridis_c()


# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Community Ordination and Environmental Gradients
# ============================================================

library(vegan)
data("mite", "mite.env")

# The mite datasets contain information on Oribatid mite communities
# sampled from a small peatland area (2.5 m × 10 m).
#
# There are linked datasets:
# ------------------------------------------------------------
# mite     : Species abundance data (35 mite species × 70 sites)
# mite.env : Environmental variables measured at the same sites
# ------------------------------------------------------------
#
# Environmental variable descriptions (mite.env):
# ------------------------------------------------------------
# SubsDens : Substrate density (g/L)
# WatrCont : Water content of the substrate (g/L)
# Substrate: Substrate type (factor with multiple levels)
# Shrub    : Shrub density (ordered factor: low → high)
# Topo     : Microtopography (Blanket vs Hummock)
# ------------------------------------------------------------

# 1. Explore and visualize interrelationships among species abundances.
#    - Examine patterns of co-occurrence.
#    - Assess whether relationships among species appear linear or nonlinear.

ggpairs(mite,
        columns = 1:5,
        progress = FALSE)

ggpairs(wisconsin(mite),
        columns = 1:5,
        progress = FALSE)


# 2. Fit a redundancy analysis (RDA) model using environmental variables of your choice.
#    - Visualize the ordination results.
#    - Examine gradients and species–environment relationships.
#    - Evaluate whether the assumptions of RDA are appropriate for these data.

df_env <- mite.env %>% 
  janitor::clean_names() %>% 
  mutate(num_shrub = as.numeric(shrub))

m_y <- mite
obj_rda_mite <- rda(m_y ~ subs_dens + watr_cont + num_shrub,
                    data = df_env)

scores(obj_rda_mite,
       display = "sites",
       scaling = 2) %>% 
  bind_cols(df_env) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%  
  ggplot(aes(x = rda1,
             y = rda2)) +
  geom_point() +
  theme_bw()

# 3. Apply alternative ordination methods.
#    - Canonical correspondence analysis (CCA; see ?cca()).
#    - Distance-based RDA (dbRDA).

obj_cca_mite <- cca(m_y ~ subs_dens + watr_cont + num_shrub,
                    data = df_env)

scores(obj_cca_mite,
       display = "sites",
       scaling = 2) %>% 
  bind_cols(df_env) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%  
  ggplot(aes(x = cca1,
             y = cca2)) +
  geom_point() +
  theme_bw()

obj_db_mite <- dbrda(m_y ~ subs_dens + watr_cont + num_shrub,
                     data = df_env,
                     distance = "bray")

scores(obj_db_mite,
       display = "sites",
       scaling = 2) %>% 
  bind_cols(df_env) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%  
  ggplot(aes(x = db_rda1,
             y = db_rda2)) +
  geom_point() +
  theme_bw()

# 4. Compare RDA, CCA, and dbRDA.
#    - Perform permutation analysis to examine the significance of predictor variables
#    - Discuss which method is most appropriate for these data and why.

anova.cca(obj_rda_mite,
          by = "margin",
          permutations = 999)

anova.cca(obj_cca_mite,
          by = "margin",
          permutations = 999)

anova.cca(obj_db_mite,
          by = "margin",
          permutations = 999)

# try variable transformation
m_mite_trans <- vegan::wisconsin(mite)

