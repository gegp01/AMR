# Requires: i = "family of antibiotics" , p = "predictor = country", and data.frame d
# libraries: rgdal, sp, sf
#p = p
#d = d
#i = i

f.plot.scores = function(x){
  i = names(d.m)[x] # variable 1
       N = as.data.frame(table(d[,p], d.m[,i]))
       names(N) = c("j", "i", "N")

  #  Likelihoods S_x
       p1 = N$N[N$i == 1]/sum(N$N[N$i == 1]) #  P(S_x | R_i)
       p2 = N$N[N$i == 0]/sum(N$N[N$i == 0]) #  P(S_x | R_noi)

  # Priors
  prior1 = sum(N$N[N$i == 1]) / sum(N$N) # P(R_i)
  prior2 = sum(N$N[N$i == 0]) / sum(N$N) # P(R_noi)
  
  scores = as.data.frame(table(d[,p], d.m[,i]))
  scores$p1 = round(p1, 3)
  scores$p2 = round(p2, 3)
  scores$score = log(((p1+0.0001)/(p2+0.0005))+log(prior1/prior2)
  #scores$score = log(((p1+0.0001)/(p2+0.0005)))

  score=scores$score
  names(score) = scores$Var1

  ncolors= 5
  x.vector=score[score>0]
  my_palette <- grDevices::colorRampPalette(c("purple", "red"))(ncolors)
  mypal=my_palette[cut(x.vector,ncolors)]
  names(mypal) = names(x.vector)

  par(mai=c(0,0,1,0))
  plot(land, col=mypal[land@data$ADMIN])
  legend("top", cex=0.7, box.col = "azure", bg="azure", horiz = T, inset = 0
       , title=paste("Resistance to", i, "(scores > 0)")
       , legend=round(c(quantile(x.vector, probs= (1:5)/5)), 2), fill = my_palette)
  }
