# GET DATA
require(ape)
require(Hmisc)
# path = "~/sandbox/AMR/"
# setwd(path)
#d = read.csv("~/sandbox/AMR/Global_AMR_enriched_FINAL_corrected.csv", sep = ";")
d = read.csv("~/sandbox/AMR/Global_AMR_enriched_FINAL_ISO3_abnamescorrected.csv", sep = ",")
d.back = d

#names(d)[names(d)=="Betalactamics"]<-"Betalactans"
#names(d)[names(d)=="Sulfonamides"]<-"Sulphonamides"


#d = d[is.na(d$year)==F,] # Eliminate data without year

antibiotics = c("Aminoglycosides","Beta.lactams","Polymyxins"
                ,"Fosfomycin","Glycopeptides","Macrolides","Oxazolidinones","Phenicols"      
                ,"Quinolones","Rifampicin","Sulphonamides","Tetracyclines","Trimethoprim")

d$amr = rowSums(d[,antibiotics])
d = d[is.na(rowSums(d[, antibiotics]))==F,] # ELIMINATE NAs IN ANTIBIOTICS


# ORGANIZE DATA: ESTIMATE MDR % PER SPECIES PER COUNTRY (CHECK THAT WE CONSIDER DATA WITH YEAR)

names(d)
d$mdr = ifelse(d$amr>=3, 1, 0)

# Split data by species
S = split.data.frame(d, d$species)

focal.taxa= c("Escherichia coli","Klebsiella pneumoniae","Staphylococcus aureus", "Enterococcus faecium", "Acinetobacter baumannii", "Pseudomonas aeruginosa")

focal.taxa = unique(d$species)

f.mdr = function(x){
  s = focal.taxa[x]
  q = aggregate(S[[s]]$mdr, list(S[[s]]$ISO.3), mean)
  names(q) =  c("ISO_3","mdr")
  q$taxa = s
  z = table(S[[s]]$ISO.3)
  sample.size = as.numeric(z)
  names(sample.size) = names(z)
  d_ = data.frame(q, sample.size)
  d_[d_$ISO_3!="",]
}

SPP = lapply(1:length(focal.taxa), f.mdr)
names(SPP) = focal.taxa

data.spp = do.call(rbind, SPP)
str(data.spp)

# NOTE: Its not possible to model mdr% ~ year, because the % is based on a set of samples.
# We can do it by a surrogate as the minimum year in the sample...



# READ PHYLOGENY OF BACTERIA SPECIES
require(ape)
V = read.tree("~/sandbox/AMR/species_level_tree.newick") 

tip_ = V$tip.label
tip_ = gsub("_", " ", tip_)

nms = d$species
nms=gsub(" ", "_", nms)
table(is.na(match(V$tip.label, nms))) # 36 specis are not in the phylogeny

# Check Which species are not in the phylogeny:
tip_[is.na(match(V$tip.label, nms))]

# 254338 samples are named as in the phylogeny, 2because they don't have at least 100 samples
table(is.na(match(nms, V$tip.label))) # And we include almost all data; 211 samples dont have a species name in the phylogeny
nms[is.na(match(nms, V$tip.label))]


# Construct phylogenetic covariance matrix
M <- vcv.phylo(V)

D = d[is.na(d$year)==F & is.na(d$country)==F,]
D$taxa = gsub(" ", "_", D$species)

# Specify model
mod0 <- mdr ~ year + ISO.3

# Construct phylogenetic covariance matrix
M <- vcv.phylo(V)

# PRIOR
prior1 <- list(
  R = list(V = 1, nu = 0.002),
  G = list(G1 = list(V = 1, nu = 0.002))
  #           , G2 = list(V = 1, nu = 0.002))  # Added G2 prior for binary response
)

# Fit model # TIME CONSUMING!!!
mod1 <- MCMCglmm(
  mod0, random = ~ taxa, family = "categorical",
  data = D, verbose = FALSE, thin = 10, burnin = 100, prior = prior1
)

summary(mod1)


########################################
# references:
# # Collect data and phylogeny
# ##########################
# # Install packages
# #install.packages("MCMCglmm")
#
# # Load library
library(MCMCglmm)
#
# # Load bird families data
data(bird.families)
#
# # Create aves data set
set.seed(123)
n_species <- length(bird.families$tip.label)
aves <- data.frame(
    especie = bird.families$tip.label,
    tamano = rnorm(n_species, 10, 1),
    sexo = sample(c("M", "F"), n_species, replace = TRUE),
    prop = rpois(n_species, lambda = 10)
    )

#Read Phylogeny
V = bird.families



# Construct phylogenetic covariance matrix
M <- vcv.phylo(V)

# Specify model with priors
mod0 <- prop ~ tamano + sexo

prior1 = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))

# Fit model
mod1 <- MCMCglmm(mod0, random = ~especie, family = "poisson", data = aves, 
                 verbose = FALSE, thin = 10, burnin = 100, prior = prior1)

# Summarize model results
summary(mod1)

###########################
## 
##  BOOLEAN - BINARY RESPONSE
##########################

# Load library
library(MCMCglmm)

# Load bird families data
data(bird.families)

# Create aves data set
set.seed(123)
n_species <- length(bird.families$tip.label)
aves <- data.frame(
  especie = bird.families$tip.label,
  tamano = rnorm(n_species, 10, 1),
  sexo = sample(c("M", "F"), n_species, replace = TRUE),
  prop = rbinom(n_species, 1, 0.5)  # Generating binary response (0 and 1)
)

# Read Phylogeny
V <- bird.families

# Construct phylogenetic covariance matrix
M <- vcv.phylo(V)

# Specify model with priors
mod0 <- prop ~ tamano + sexo


prior1 <- list(
  R = list(V = 1, nu = 0.002),
  G = list(G1 = list(V = 1, nu = 0.002))
#           , G2 = list(V = 1, nu = 0.002))  # Added G2 prior for binary response
)

# Fit model
mod1 <- MCMCglmm(
  mod0, random = ~ especie, family = "categorical",
  data = aves, verbose = FALSE, thin = 10, burnin = 100, prior = prior1
)

# Summarize model results
summary(mod1)


#################################################################
# MULTIPLE RESPONSE MODEL
# Load bird families data
data(bird.families)

# Create aves data set
set.seed(123)
n_species <- length(bird.families$tip.label)
n_rows <- 100
n_cols <- 13
aves <- data.frame(
  especie = rep(bird.families$tip.label, each = n_rows),
  tamano = rnorm(n_species * n_rows, 10, 1),
  sexo = rep(sample(c("M", "F"), n_species, replace = TRUE), each = n_rows),
  gen = matrix(sample(c(0, 1), n_species * n_rows * n_cols, replace = TRUE), nrow = n_species * n_rows)
)


especie = rep(bird.families$tip.label, each = n_rows)
tamano = rnorm(n_species * n_rows, 10, 1)
sexo = rep(sample(c("M", "F"), n_species, replace = TRUE), each = n_rows)
gen = matrix(sample(c(0, 1), n_species * n_rows * n_cols, replace = TRUE), nrow = n_species * n_rows)


# Construct phylogenetic covariance matrix
M <- vcv.phylo(bird.families)

# Specify model with priors
mod0 <- gen ~ tamano + sexo
prior1 = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))

# Fit model
mod1 <- MCMCglmm(mod0, random = ~especie, family = "multinomial", data = aves,
                 verbose = FALSE, thin = 10, burnin = 100, prior = prior1)

# Summarize model results
summary(mod1)

