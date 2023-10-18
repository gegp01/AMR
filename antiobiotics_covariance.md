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

#### Covariance model MCMCglmm
Monte Carlo Simulations require RAM, hence analysis of large dataset may crash the computer.
A work around could be to analyse part of the data to train a model; and then the model can be validated by testing its predictions on the other part of the dataset.

##### Select a random sample of the dataet to train a model.
~~~~
set.seed(193839)
s = as.integer(nrow(D)*0.5)
D.sub.tr = D[sample(1:nrow(D), s),]

set.seed(193839)
D_n0 = D[rowSums(D[,antibiotics])!=0,]
s = as.integer(nrow(D_n0)*0.5)
D.sub.tr.n0 = D_n0[sample(1:nrow(D_n0), s),]

~~~~


#### Set priors
As prior belief for the correlations among resistance to antibiotic families, we can 
use a matrix V of 13 x 13 with an expected self-correlation of 1, and no correlation between resistance to different antibiotic families. In other words, the elements in the diagonal of the matrix V are 1, and the elements off-diagonal are 0.

This priors will be challenged against our hypothesis and data.

~~~
# For one factor
usP_1<-list(V = diag(13), nu = 2)

# For two factors we can split the variance 
usP_2<-list(V = diag(13)/5, nu = 2)

# For three factors
usP<-list(V = diag(13)/3, nu = 2)
idhP<-list(V = diag(13)*2/3, nu = 2)

~~~
##### Run de model
May take some time...
~~~
#LIBRARY
require(MCMCglmm)

# PRIOR
priorX <- list(R = identity, G = list(G1 = identity, G2 = identity))

# RUN THE MODEL: CATEGORICAL MULTIVARIATE RESPONSE WITH A COVARIATE MATRIX OF ANTIBIOTICS

modelX_rnd <-MCMCglmm(cbind(Aminoglycosides,Beta.lactams,Polymyxins,Fosfomycin
                         ,Glycopeptides,Macrolides,Oxazolidinones,Phenicols
                         ,Quinolones,Rifampicin,Sulphonamides
                         ,Tetracyclines,Trimethoprim) ~ 1 #+ country - 1
                        , family=rep("categorical", 13)
                   , random = ~ us(trait):taxa + us(trait):country
                   , rcov = ~ us(trait):units
                   , prior=priorX
                   , data= D.sub.tr.n0)
~~~

##### Call the covariance matrix and plot the bipartite network.

~~~
modelX4$VCV

plot(network)
~~~

