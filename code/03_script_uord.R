#' DESCRIPTION:
#' Script for Unconstrained Ordination

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               GGally,
               vegan)

## PCA ####
## use iris data
df_iris <- iris %>% 
  as_tibble() %>% 
  janitor::clean_names()

df_iris %>% 
  ggpairs(progress = FALSE,
          columns = c("sepal_length",
                      "sepal_width",
                      "petal_length",
                      "petal_width"),
          aes(color = species,
              alpha = 0.5)
          ) +
  theme_bw()

## starts_with() selects columns based on key
df_petal <- df_iris %>% 
  select(starts_with("petal_"))

## PCA - prcomp() function
obj_pca <- prcomp(x = df_petal,
                  center = TRUE,
                  scale = TRUE)

## get more information
summary(obj_pca)

## get PC axes values
df_pca <- df_iris %>% 
  bind_cols(obj_pca$x)

## draw a figure comparing PC1 values between species
df_pca %>% 
  ggplot(aes(x = species,
             y = PC1)) +
  geom_boxplot()


## NMDS ####

# call dune data
data(dune)

# visual
dune %>% 
  as_tibble() %>% 
  select(1:3) %>% 
  ggpairs() + 
  theme_bw()

# calculate distance between units
# - column is species
# - row is site (or similar)
m_bray <- vegdist(dune, method = "bray")

obj_nmds <- metaMDS(comm = m_bray,
                    k = 2)

# metaMDS(comm = dune,
#         distance = "bray",
#         k = 2)

# visualize NMDS
data(dune.env)

df_nmds <- dune.env %>% 
  as_tibble() %>% 
  bind_cols(obj_nmds$points) %>% 
  janitor::clean_names()

df_nmds %>% 
  ggplot(aes(x = mds1,
             y = mds2,
             color= use)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95,
               linetype = 2) +
  theme_bw() +
  labs(x = "NMDS1",
       y = "NMDS2",
       color = "Land use intensity")

# permanova
adonis2(m_bray ~ use,
        data = df_nmds)

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: PCA using the iris dataset
# ============================================================

# In this exercise, you will perform a Principal Component
# Analysis (PCA) using all morphological measurements in the
# iris dataset and visualize multivariate trait patterns
# among species.

# 1. Using all four morphological variables
#    (Sepal.Length, Sepal.Width, Petal.Length, Petal.Width),
#    perform a PCA.

df_ps <- df_iris %>% 
  select(-species)

obj_pca4 <- prcomp(x = df_ps,
                   scale = TRUE,
                   center = TRUE)

summary(obj_pca4)

df_pca4 <- df_iris %>% 
  bind_cols(obj_pca4$x)

# 2. Visualize flower morphology in PC axes whose cumulative
#    contribution exceeds 90%; color points by species.

df_pca4 %>% 
  ggplot(aes(x = PC1,
             y = PC2,
             color = species)) +
  geom_point()

# 3. Which morphological traits contribute most strongly to
#    the first and second principal components? How?

# visualize loadings of orginal variables
loadings <- obj_pca4$rotation[, 1:2]
eig <- obj_pca4$sdev^2
loadings_scaled <- sweep(loadings, 2, sqrt(eig[1:2]), "*")

df_pca4 %>% 
  ggplot(aes(x = PC1,
             y = PC2,
             color = species)) +
  geom_point() +
  geom_segment(data = loadings_scaled,
               aes(x = 0,
                   xend = PC1,
                   y = 0,
                   yend = PC2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm"))
               )

# ============================================================
# EXERCISE: NMDS using the BCI dataset
# ============================================================

# In this exercise, you will perform a Non-metric Multidimensional
# Scaling (NMDS) using the BCI tree community dataset and explore
# patterns in species composition among sites.

data("BCI", "BCI.env")

# 1. Using the BCI dataset, calculate a dissimilarity matrix
#    (e.g., Bray-Curtis) and perform NMDS.

m_bray_bci <- vegdist(BCI, method = "bray")
obj_nmds_bci <- metaMDS(comm = m_bray_bci,
                        k = 2)

# 2. Visualize sites in NMDS space.
#    - How are sites positioned relative to each other?
#    - Color or shape points by environmental groups or site
#      characteristics of your choice.

df_nmds_bci <- BCI.env %>% 
  as_tibble() %>% 
  bind_cols(obj_nmds_bci$points) %>% 
  janitor::clean_names()

df_nmds_bci %>% 
  ggplot(aes(x = mds1,
             y = mds2,
             color = habitat)) +
  geom_point(size = 3)

# 3. Perform PERMANOVA to examine if communities are grouped
#    by the environmental variable you selected.
