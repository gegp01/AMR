# Objectives:
# Estimate the posterior distribution of the parameters that explain statistical differences between clases.
# The clasess to investigate are:
# (1) Differences between species within the same country or region.
# (2) Differences between countries, within each species.
# These functions herein will select a balanced sample of species within the region.
# Conditions for the data subset:
# A balanced sample of data will have the same number of observations for all species.
# The data subset must contain similar ammount of zeros and ones in the subset.

path = "PATH TO THE DATA"
setwd(path)

# 1. Read data
D= read.csv("Global_AMR_Final_Database_2_corrected_coordinates 1.csv")

# 2. Select Regions
require(rgdal)
land = readOGR("https://gegp01.github.io/ServSoc/countries.geojson")

# 2.1 Read species phylogeny
require(ape)
V = read.tree("https://gegp01.github.io/AMR/")



