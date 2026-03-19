#' DESCRIPTION:
#' Script for time-series

# in-class ----------------------------------------------------------------

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
