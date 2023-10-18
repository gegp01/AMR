### Correlations of resistance between antibiotic families.
Gabriel E García Peña 16-oct-2023

#### Data
~~~
# LIBRARIES
require(ape)
require(Hmisc)
require(AMR)
require(dplyr)

# DATA
d = read.csv("~/sandbox/AMR/Global_AMR_enriched_FINAL_ISO3_abnamescorrected.csv", sep = ",")

# CLASSIFY SPECIES AS GRAM + OR GRAM -
d_mo<-d%>%mutate(bacteria=as.mo(species))#data es el nombre de la base de datos
d_gram<-d_mo%>%mutate(gram=mo_gramstain(bacteria))

d = d_gram

antibiotics = c("Aminoglycosides","Beta.lactams","Polymyxins"
                ,"Fosfomycin","Glycopeptides","Macrolides","Oxazolidinones","Phenicols"      
                ,"Quinolones","Rifampicin","Sulphonamides","Tetracyclines","Trimethoprim")

# ESTIMATE AMR
d$amr = rowSums(d[,antibiotics])
d = d[is.na(rowSums(d[, antibiotics]))==F,] # ELIMINATE NAs IN ANTIBIOTICS

# CLASSIFY MDR AS AMR > 3
d$mdr = ifelse(d$amr>=3, 1, 0)
~~~

#### Curate Data
Location and year may be missing for some data. Also, there are some misstakes that must be fixed.
As working data we can consider the whole dataset, or the dataset 

~~~
# SELECT SPECIES WITH AT LEAST 100 SAMPLES GLOBALY. THIS IS OPTIONAL.
# TO PERFORM THE ANALYSIS ON THE WHOLE DATASET JUST CHANGE d_ = d_100 TO d_ = d  
X = table(d$species)
taxa = names(X[X>=100])

d.100 = d[is.na(match(d$species, taxa))==F,]
d_ = d.100

# CLEAN GEOGRAPHIC DATA. CURATE COUNTRIES AND ISO3
D = d_[is.na(d_$year)==F & is.na(d_$country)==F,]
D$taxa = gsub(" ", "_", D$species)

D$ISO.3[D$country=="Morocco"]<-"MAR"
D$country[D$country=="Morocco"]<-"Morocco"
D$country[grep("UK", D$country)]<-"United Kingdom"
D$country[grep("Russia", D$country)]<-"Russia"
D$country[grep("Unknown", D$country)]<-NA

# EXCLUDE DATA WITHOUT COUNTRY
D = D[is.na(D$country)==F,]

# FILL MISSING REGIONS
D$continent[grep("Argentina", D$country)]<-"South America"
D$continent[grep("Brazil", D$country)]<-"South America"
D$continent[grep("Germany", D$country)]<-"Europe"
D$continent[grep("United Kingdom", D$country)]<-"Europe"
D$continent[grep("Italy", D$country)]<-"Europe"
D$continent[grep("Russia", D$country)]<-"Asia"
D$continent[grep("Cameroon", D$country)]<-"Africa"
D$continent[grep("Belgium", D$country)]<-"Europe"
D$continent[grep("Trinidad and Tobago", D$country)]<-"South America"
D$continent[grep("Burkina Faso", D$country)]<-"Africa"
~~~


#### Phylogenetic signal
 MDR
 AMR
 AND THE WHOLE matrix {0,1}

#### Model MCMCglmm on the correlations between resistance to antibiotic families.
The model presented here aims at capturing the posterior of the correlations between resistance to one of the 13 antibiotic families and resistance to any of the other antibiotic families. 

##### Select a random sample of the dataet to train a model.
Monte Carlo Simulations require RAM, hence analysis of large dataset may crash the computer. A work around could be to analyse part of the data to train a model; and then the model can be validated by testing its predictions on the other part of the dataset.

~~~~
set.seed(193839)
s = as.integer(nrow(D)*0.5)
D.sub.tr = D[sample(1:nrow(D), s),]

set.seed(193839)
D_n0 = D[rowSums(D[,antibiotics])!=0,]
s = as.integer(nrow(D_n0)*0.5)
D.sub.tr.n0 = D_n0[sample(1:nrow(D_n0), s),]

~~~~


##### Set priors
As prior belief for the correlations among resistance to antibiotic families, we can use a matrix of 13 x 13 with an expected self-correlation of 1, and no correlation between the resistance to different antibiotic families. In other words, the elements in the diagonal of the matrix are 1, and the elements off-diagonal are 0. 

This prior will be tested against the data; and we can modifify the prior values of self-correlation as in the following examples.
Furthermore, we could change the elements of diagonal if there is some hypothesis about it.

~~~
# Prior matrix considering self-correlation of 1.
identity<-list(V = diag(13), nu = 2)

# Prior matrix considering self-correlation of 0.33.
var_33<-list(V = diag(13)/3, nu = 2)

# Prior matrix considering self-correlation of 0.66.
var_66<-list(V = diag(13)*2/3, nu = 2)

~~~
##### Run de model
May take some time...
~~~
#LIBRARY
require(MCMCglmm)

# SET PRIOR FOR THE CORRELATION MATRIX (rcov = corg(traits):units)
prior1 <- list(R = identity)

# RUN THE MODEL: MULTIVARIATE CATEGORICAL RESPONSE WITH A COVARIATE MATRIX OF ANTIBIOTICS

require(MCMCglmm)
modelX_corg <-MCMCglmm(cbind(Aminoglycosides,Beta.lactams,Polymyxins,Fosfomycin
                            ,Glycopeptides,Macrolides,Oxazolidinones,Phenicols
                            ,Quinolones,Rifampicin,Sulphonamides
                            ,Tetracyclines,Trimethoprim) ~ 1 #+ country - 1
                      , family=rep("categorical", 13)
                      , rcov = ~ corg(trait):units # FIXING THE COVARIANCE MATRIX TO BE A CORRELATION MATRIX
                      , prior=prior1
                      , data= D.sub.tr.n0)

~~~

##### Call the correlation matrix and plot the bipartite network.

~~~
modelX4$VCV

plot(network)
~~~

