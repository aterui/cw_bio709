#' DESCRIPTION:
#' Script for SEM

# in-class ----------------------------------------------------------------

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Path Analysis and Covariance Visualization
# ============================================================

library(piecewiseSEM)
data("keeley")

# The "keeley" dataset contains fire-related vegetation data
# collected from shrublands in California.
#
# ------------------------------------------------------------
# Column descriptions:
# elev  : Elevation of the site
# abiotic : Abiotic environmental factors (integrated variable indicating the suitability of the habitat)
# age   : Stand age
# hetero  : Environmental heterogeneity
# firesev: Fire severity
# cover : Vegetation cover
# rich  : Plant species richness
# ------------------------------------------------------------
#
# In this exercise, you will explore relationships among variables
# using covariance and path analysis. You will replicate a published
# path model and propose an alternative.

# 1. For the variables depicted in Figure 22.1, draw a figure
#    showing the covariance between variables.

# 2. Following Figure 22.1, develop a path model using the
#    same variables and relationships. Examine if this model
#    captures the data structure using a Chi-Square test.

# 3. Develop an alternative path model that you consider more
#    appropriate based on theory or observed data patterns.

# 4. Compare the performance of the published model (Figure 22.1)
#    and your alternative model.
#    - Consider fit indices, path coefficients, and interpretability.
