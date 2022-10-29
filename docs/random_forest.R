# BOOTSTRAP TREE AND DATA ANALYSIS
f0 = function(x) {source("https://gegp01.github.io/AMR/bootstrap.R")}
G = lapply(1:100, f0)


# Collect bootstrap for AUSTRALIA
f.data = function(x) {
  AMR = names(G[[1]][[1]][[1]][[1]][,1:13])
  com = G[[x]][[1]][[1]][[1]][,AMR]
  row.names(com) = gsub(" ", "_", G[[1]][[1]][[1]][[1]]$species)
  rowSums(com)
}

AMR.AUS = lapply(1:length(G), f.data)


# Collect bootstrap for CHINA
f.data = function(x) {
  AMR = names(G[[1]][[1]][[1]][[2]][,1:13])
  com = G[[x]][[1]][[1]][[2]][,AMR]
  row.names(com) = gsub(" ", "_", G[[1]][[1]][[1]][[2]]$species )
  rowSums(com)
}


AMR.CHINA = lapply(1:length(G), f.data)


# Collect bootstrap for UK
f.data = function(x) {
  AMR = names(G[[1]][[1]][[1]][[3]][,1:13])
  com = G[[x]][[1]][[1]][[3]][,AMR]
  row.names(com) = gsub(" ", "_", G[[1]][[1]][[1]][[3]]$species)
  rowSums(com)
}

AMR.UK = lapply(1:length(G), f.data)


# Collect bootstrap for USA
f.data = function(x) {
  AMR = names(G[[1]][[1]][[1]][[4]][,1:13])
  com = G[[x]][[1]][[1]][[4]][,AMR]
  row.names(com) = gsub(" ", "_", G[[1]][[1]][[1]][[4]]$species)
  rowSums(com)
}

AMR.USA = lapply(1:length(G), f.data)


# Box plot comparison between species

pdf("species_AMR.pdf")
par(mai=c(3, 1, 0.5, 0.5), mfrow=c(1,1), font.main=1)
  boxplot(do.call(rbind, AMR.AUS), las= 2, cex= 0.5, main="Australia", ylab="AMR")
  boxplot(do.call(rbind, AMR.CHINA), las= 2, cex= 0.5, main= "China", ylab="AMR")
  boxplot(do.call(rbind, AMR.UK), las= 2, cex= 0.5, main="United Kingdom", ylab="AMR")
  boxplot(do.call(rbind, AMR.USA), las= 2, cex= 0.5, main= "USA", ylab="AMR")
dev.off()

