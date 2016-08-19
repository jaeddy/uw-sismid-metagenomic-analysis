
## 1. Load the data
load('STAT.RData')
library(phyloseq)


head(sample_data(phy))
levels(sample_data(phy)$Treatment)

head(phenotypes)
phy
sample_names(phy)[1:5]

## The sample names are in different formats
## phenotypes are per mouse
## phy data are per mouse-location

## 2. Merge phenotypes
pheno1 = phenotypes
pheno2 = phenotypes

rownames(pheno1) = paste('cecal', rownames(pheno1), sep="_")
rownames(pheno2) = paste('fecal', rownames(pheno2), sep="_")
pheno = rbind(pheno1, pheno2)

phy2 = merge_phyloseq(phy, sample_data(pheno))

phy
phy2
head(sample_data(phy2))

## cleanup
rm(pheno, pheno1, pheno2, phy, phenotypes)

## Estimate richness
?estimate_richness

alpha = estimate_richness(phy2, measures=c("Observed", "Chao1", "Shannon"))
head(alpha)

plot(ecdf(sample_sums(phy2)))


### Let's only keep samples that have at least 3500 observations/reads
phy3 = subset_samples(phy2, sample_sums(phy2)>3500)
alpha = estimate_richness(phy3, measures=c("Observed", "Chao1", "Shannon"))

alpha.rare = matrix(0, ncol=ncol(alpha), nrow=nrow(alpha))
for(i in 1:20){
  alpha.rare= alpha.rare + 
    estimate_richness(rarefy_even_depth(phy3, sample.size = 3500,
                                        verbose=F, trimOTUs = F), 
                      measures=c("Observed", "Chao1", "Shannon"))
}
alpha.rare = alpha.rare/20

## Compare the rarefied vs non-rarefied alpha diversity
plot(alpha$Observed, alpha.rare$Observed)
plot(alpha$Chao1, alpha.rare$Chao1)
plot(alpha$Shannon, alpha.rare$Shannon)

## Compute order level abundances
order.phy = tax_glom(phy3, taxrank = "Order")

## Normalize the data

# relative abundance
order.rel = transform_sample_counts(order.phy, function(x) x/sum(x))

# clr
order.clr = transform_sample_counts(order.phy, function(x){y=log(x+1); y/sum(y)})

# DESeq2
library(DESeq2)
countData = round(as(otu_table(order.phy), "matrix"), digits = 0)
countData = countData + 1L
dds <- DESeqDataSetFromMatrix(countData, sample_data(order.phy), design = ~1)

order.deseq = merge_phyloseq(otu_table(counts(estimateSizeFactors(dds), 
                                              normalized=T), taxa_are_rows = T), 
                             sample_data(order.phy), 
                             phy_tree(order.phy), 
                             tax_table(order.phy))

## cleanup
rm(dds, countData)

# Compare location
Location = sample_data(order.rel)$Location
rel.p = apply(otu_table(order.rel),1, function(x) wilcox.test(c(x)~Location)$p.value)
clr.p = apply(otu_table(order.clr),1, function(x) t.test(c(x)~Location)$p.value)
deseq.p = apply(otu_table(order.deseq),1, function(x) t.test(log10(c(x))~Location)$p.value)

univ.location.res = data.frame(taxa = tax_table(order.clr)[,"Order"], rel.p, clr.p, deseq.p)

# do FDR
univ.location.res$rel.fdr = p.adjust(univ.location.res$rel.p, method="fdr")
univ.location.res$clr.fdr = p.adjust(univ.location.res$clr.p, method="fdr")
univ.location.res$deseq.fdr = p.adjust(univ.location.res$deseq.p, method="fdr")

colSums(univ.location.res<0.05)
univ.location.res

## Plot heatmaps
library(ggplot2)
plot_heatmap(order.rel, taxa.label = "Order", sample.order = "Location")+
  facet_grid(.~Location, scales = "free_x")
plot_heatmap(order.clr, taxa.label = "Order", sample.order = "Location")+
  facet_grid(.~Location, scales = "free_x")
plot_heatmap(order.deseq, taxa.label = "Order", sample.order = "Location")+
  facet_grid(.~Location, scales = "free_x")

## Plot barplots of all data
plot_bar(order.rel, fill="Location") + 
  facet_wrap(facets = "Order", ncol=7, nrow=3, scales = "free_y")
plot_bar(order.clr, fill="Location") + 
  facet_wrap(facets = "Order", ncol=7, nrow=3, scales = "free_y")
plot_bar(order.deseq, fill="Location") + 
  facet_wrap(facets = "Order", ncol=7, nrow=3, scales = "free_y")

## Plot box and whiskers charts
ggplot(psmelt(order.rel), aes(x=Location, y=Abundance)) + 
  geom_boxplot(aes(fill=Location)) + geom_jitter() + 
  facet_wrap(facets="Order", ncol=7, nrow=3, scales="free_y")

ggplot(psmelt(order.clr), aes(x=Location, y=Abundance)) + 
  geom_boxplot(aes(fill=Location)) + geom_jitter() + 
  facet_wrap(facets="Order", ncol=7, nrow=3, scales="free_y")

ggplot(psmelt(order.deseq), aes(x=Location, y=Abundance)) + 
  geom_boxplot(aes(fill=Location)) + geom_jitter() + 
  facet_wrap(facets="Order", ncol=7, nrow=3, scales="free_y") + scale_y_log10()

## Fecal samples vs Treatment
fecal.rel = subset_samples(order.rel, Location = "fecal")
fecal.rel.subs = subset_taxa(fecal.rel, rowMeans(otu_table(fecal.rel))> 0.001)
Treatment = sample_data(fecal.rel.subs)$Treatment

fecal.res = data.frame(taxa = tax_table(fecal.rel.subs)[,"Order"], 
                               fecal.p = apply(otu_table(fecal.rel.subs),1, 
                                               function(x) kruskal.test(c(x)~Treatment)$p.value))
fecal.res$fecal.fdr = p.adjust(fecal.res$fecal.p, method="fdr")

ggplot(psmelt(fecal.rel.subs), aes(x=Treatment, y=Abundance)) + 
  geom_boxplot(aes(fill=Treatment)) + geom_jitter() + 
  facet_wrap(facets="Order", ncol=4, nrow=2, scales="free_y")


## Fecal samples vs BMD
fecal.deseq = subset_samples(order.deseq, Location = "fecal")

pheno.res = data.frame(taxa = tax_table(fecal.deseq)[,"Order"], 
                       BMD.p = apply(otu_table(fecal.deseq),1, 
                                       function(x) cor.test(log(c(x)), sample_data(fecal.deseq)$BMD)$p.value))
pheno.res$BMD.fdr = p.adjust(pheno.res$BMD.p, method="fdr")

pheno.res$pFAT.p = apply(otu_table(fecal.deseq),1, 
                         function(x) cor.test(log(c(x)), sample_data(fecal.deseq)$pFat)$p.value)

pheno.res$pFAT.fdr = p.adjust(pheno.res$pFAT.p, method="fdr")

gFDR = with(pheno.res, 
            matrix(p.adjust(c(BMD.p, pFAT.p), method="fdr"), 
                   ncol=2, nrow=nrow(pheno.res), byrow=F))

gFDR

## Plot phenotypic variables vs. "the most significant" result
ggplot(subset(psmelt(fecal.deseq), Order=="Enterobacteriales"), aes(x=BMD, y=Abundance)) +
  geom_point() + scale_y_log10()
ggplot(subset(psmelt(fecal.deseq), Order=="Enterobacteriales"), aes(x=pFat, y=Abundance)) + 
  geom_point() + scale_y_log10()


## Plot phenotypic variables vs. Treatment
ggplot(sample_data(fecal.deseq), aes(x=Treatment, y=BMD)) + 
  geom_boxplot(aes(fill=Treatment))+geom_jitter()
ggplot(sample_data(fecal.deseq), aes(x=Treatment, y=pFat)) + 
  geom_boxplot(aes(fill=Treatment))+geom_jitter()




