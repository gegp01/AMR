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
V = read.tree("https://gegp01.github.io/AMR/SpeciesLevelTree.newick")
nms_tree = V$tip.label


# DATA
d.back = D # Make a back up of the original data.
D = D[is.na(d$year)==F,] # Eliminate data without year

# Change name of antibiotic family
names(d.sample)[names(d.sample)=="amino"]<-"aminoglycoside"

# Names of the antibiotics
antibiotics = c("aminoglycoside", "betalactamics", "colistin", "fosfomycin", "glycopeptide", "macrolide"
                , "oxazolidinone", "phenicol", "quinolone", "rifampicin", "sulfo0mide", "tetracycline", "trimethoprim")

# SUBSET MAKER

# 1. DETERMINE MINIMUM NUMBER OF SAMPLES PER SPECIES
min.sample = 600

# DEFINE TIME WINDOW (YEAR)
year.start = 2010
year.end = 2022

############################################# INLCUDE THIS IN A FUNCTION
# SUBSET DATASET BASED ON SAMPLE SIZE
X = table(d$species)
taxa = names(X[X>=min.sample) # select taxa

d.sample = d[is.na(match(d$species, taxa))==F,]


# SUBSET OF CURRENT TIMES (las 10 years)
d.sample.time = d.sample[d.600$year>=year.start,]

# WORKING DATASET
D = d.sample.time
##############################################
##############################################
t.min = year.start
s.min = min.sample
               
f1 = function(t.min, s.min) {
  X = table(d$species)
  taxa = names(X[X>=s.min) # select taxa
  d.sample = d[is.na(match(d$species, taxa))==F,]
                 
  # SUBSET OF CURRENT TIMES (las 10 years)
  d.sample.time = d.sample[d.sample$year>=t.min,]
                 
  # WORKING DATASET
  D = d.sample.time
  D
  }


               
               

