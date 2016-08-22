library(phyloseq)
load("STAT.RData")

## Compute JSD distance matrix
jsd.dist = phyloseq::distance(phy, "jsd")

## Recall the PCoA results
jsd.pco = dudi.pco(cailliez(jsd.dist), scannf=F, nf=2)

s.class(jsd.pco$li, sample_data(phy)$Location)
s.class(jsd.pco$li, sample_data(phy)$Treatment)
s.class(jsd.pco$li, with(sample_data(phy), 
                         interaction(Location,Treatment)))

##
?adonis
samp.data = data.frame(sample_data(phy))

adonis(jsd.dist ~ Location+Treatment, data=samp.data, 
       permutations=9999)

adonis(jsd.dist ~ Treatment, data=samp.data, 
       strata = samp.data$Location)


jsd.wca = wca(jsd.pco, samp.data$Location, scannf=F, nf=2)

s.class(jsd.wca$li, samp.data$Location)
s.class(jsd.wca$li, samp.data$Treatment)

## Post-hoc tests
### Function for post-hoc test on two levels
## with stratification
Treatment = samp.data$Treatment
dist.mat=as.matrix(jsd.dist)
phoc.test = function(x){
  spair = (samp.data$Treatment %in% x)
  phoc = adonis(as.dist(dist.mat[spair,spair]) ~ Treatment[spair], 
                data=samp.data, strata = samp.data$Location[spair])
  c(x, phoc$aov.tab[,6][1])
}
## Run the function above for all pairs of Treatment levels
combn(levels(Treatment),2, phoc.test )

### Function for post-hoc test on two levels
## without stratification
phoc.test2 = function(x){
  spair = (samp.data$Treatment %in% x)
  phoc = adonis(as.dist(dist.mat[spair,spair]) ~ Treatment[spair], 
                data=samp.data)
  c(x, phoc$aov.tab[,6][1])
}
combn(levels(Treatment), 2, phoc.test2 )


###
##
###
sample_data(phy)$Abx = factor(Treatment != 'C')
fecal.phy = subset_samples(phy, Location == 'fecal')
fecal.dist = phyloseq::distance(fecal.phy, method="jsd")

adonis(fecal.dist~sample_data(fecal.phy)$Abx)


simulate.ps = function(){
  n1 = 10; x1 = rnorm(n1, mean=0, sd = 10)
  n2 = 36; x2 = rnorm(n2, mean=0, sd = 1)
  x=c(x1,x2)
  f=factor(c(rep(1,n1), rep(2,n2)))
  dm = dist(x)
  c(adonis.p = adonis(dm~f)$aov.tab[1,6],
  tev.p = t.test(x1, x2, var.equal = T)$p.value,
  tuv.p = t.test(x1, x2, var.equal = F)$p.value)
}

set.seed(20160726)
simulate.ps()
sim = t(replicate(1000, simulate.ps()))
colSums(sim<0.05)/1000


###  Tw2 statistic
WT = function(dm, f){
  if(nlevels(f) != 2) return(NULL)
  lev = levels(f)
  ns = table(f)
  N = sum(ns)
  dd = as.matrix(dm)
  dd[upper.tri(dd)]=0 ##
  SST = sum(dd^2)/N
  SSW1 = sum(dd[f==lev[1],f==lev[1]]^2)/ns[1] 
  SSW2 = sum(dd[f==lev[2],f==lev[2]]^2)/ns[2]
  SSW = SSW1 + SSW2
  
  s1 = SSW1/(ns[1]-1)
  s2 = SSW2/(ns[2]-1)
  if(SST < SSW) 
    t.stat = 0
  else
    t.stat = sqrt((ns[1]+ns[2])/(ns[1]*ns[2]))*sqrt(SST-SSW)/sqrt(s1/ns[1] + s2/ns[2])
  t.stat
}

### Permutation test using Tw2

WT.test = function(dm, f, nrep=999){
  stats = c(WT(dm, f), replicate(nrep, WT(dm, f[sample(length(f))])))
  p.value = sum(stats>=stats[1])/(nrep+1)
  t.stat = stats[1]
  list(p.value = p.value, t.stat = t.stat, nrep=nrep)  
}

WT.test(fecal.dist, sample_data(fecal.phy)$Abx)
