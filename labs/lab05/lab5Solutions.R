load('STAT.RData')
library(phyloseq)
## subset to fecal samples only
fecal.phy = subset_samples(phy, Location=='fecal')

rownames(phenotypes) = paste("fecal", rownames(phenotypes), sep="_")
fecal.phy = merge_phyloseq(fecal.phy, 
                            sample_data(phenotypes))

## compute JSD distance
jsd.fecal = phyloseq::distance(fecal.phy, method="jsd")

# compute inter-/intra-group distances for Treatment
library(vegan)
md = meandist(jsd.fecal, sample_data(fecal.phy)$Treatment)
md

## Plot the intra-group diversity
barplot(diag(md))

## PCoA
library(ade4)
dudi.pco(jsd.fecal)

## Compute PCoA
fecal.pco = dudi.pco(cailliez(jsd.fecal), 
                     scannf=F, nf=2)

## Calculate % inertia explained
fecal.pco$eig[1:2]/sum(fecal.pco$eig)


## plot communities versus factor variables
s.class(fecal.pco$li, sample_data(fecal.phy2)$Treatment)

# Plot additional phenotypes on PCoA
s.value(fecal.pco$li, sample_data(fecal.phy2)$GIP, sub="GIP")
s.value(fecal.pco$li, sample_data(fecal.phy2)$Insulin, sub="Insulin")
s.value(fecal.pco$li, sample_data(fecal.phy2)$Leptin, sub="Leptin")
s.value(fecal.pco$li, sample_data(fecal.phy2)$IGF1, sub="IGF1")
s.value(fecal.pco$li, sample_data(fecal.phy2)$BMD, sub="BMD")
s.value(fecal.pco$li, sample_data(fecal.phy2)$mFat, sub="mFat")
s.value(fecal.pco$li, sample_data(fecal.phy2)$pFat, sub="pFat")

cor.mat = with(sample_data(fecal.phy2),
               cor(cbind(PCo=fecal.pco$li,
                         GIP, Insulin, Leptin,
                         IGF1, BMD, mFat, pFat)))

round(cor.mat[,1:2], digits=2)

## plot correlations of PCo with phenotypes
library(ggplot2)
library(reshape)
mcor = melt(as.matrix(cor.mat[,1:2]))
mcor$X1 = factor(mcor$X1, levels=colnames(cor.mat))
mcor$X2 = factor(mcor$X2, levels=colnames(cor.mat))
qplot(x=X1, y=X2, 
      data=mcor, 
      fill=value, geom="tile")+
  scale_fill_gradientn(limits = c(-1,1), 
                       colours=c("red", "white", "green")) + 
  geom_text(aes(label=round(value, digits=2)))


## compute JSD distance for the entire dataset
jsd = phyloseq::distance(phy, method="jsd")

## explore the eigenvalue plot
dudi.pco(jsd)

## compute PCoA
jsd.pco = dudi.pco(cailliez(jsd), scannf=F, nf=2)

## plot communities versus factor variables
s.class(jsd.pco$li, 
        sample_data(phy)$Location)
s.class(jsd.pco$li, 
        sample_data(phy)$Treatment)
s.class(jsd.pco$li, 
        with(sample_data(phy), 
             interaction(Treatment, Location)))

jsd.wca = wca(jsd.pco, 
              sample_data(phy)$Location)

s.class(jsd.wca$li, sample_data(phy)$Location)
s.class(jsd.wca$li, sample_data(phy)$Treatment)


