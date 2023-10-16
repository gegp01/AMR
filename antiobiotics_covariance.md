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


