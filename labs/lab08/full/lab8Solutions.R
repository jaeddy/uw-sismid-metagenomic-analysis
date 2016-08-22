library(phyloseq)
load('STAT.picrust.RData')

genes

head(pathways)

library(DESeq2)
countData = round(as(otu_table(genes), "matrix"), digits = 0)
countData = countData + 1L

dds <- DESeqDataSetFromMatrix(countData, sample_data(genes), design = ~1)
genes.deseq = merge_phyloseq(log(otu_table(counts(estimateSizeFactors(dds), 
                                              normalized=T), taxa_are_rows = T)), 
                             sample_data(genes),
                             tax_table(genes))

###
# PCA
###
genes.eucl = phyloseq::distance(genes.deseq, "euclidean")

genes.pca = ordinate(genes.deseq, "PCoA", distance=genes.eucl)
plot_ordination(genes.deseq, ordination = genes.pca, color = "Location", shape="Treatment")

# Multivariate analysis
adonis(genes.eucl ~ sample_data(genes.deseq)$Location)
adonis(genes.eucl ~ sample_data(genes.deseq)$Treatment)


###
# Calculate univatiate tests
###
Location = sample_data(genes.deseq)$Location
p.values = apply(otu_table(genes.deseq), 1, function(x) t.test(c(x) ~ Location)$p.value)
q.values = p.adjust(p.values, method="fdr")

###
## Pathway enrichment
###
p = "Human Diseases"
ctab = table(names(p.values) %in% pathways$KO[pathways$Level1 == p], p.values<0.05)
fisher.test(ctab)

test.pathway.level1 = function(p){
  ctab = table(names(p.values) %in% pathways$KO[pathways$Level1 == p], p.values<0.05)
  fisher.test(ctab)$p.value
}
level1.pvalues = sapply(unique(pathways$Level1), test.pathway.level1)


## Metabolism is significant; let's look at Level2 within Metabolism
test.metabolism.pathways = function(p){
  ctab = table(names(p.values) %in% 
                 pathways$KO[pathways$Level1 == 'Metabolism' & pathways$Level2 == p], 
               p.values<0.05)
  fisher.test(ctab)$p.value
}
metabolism.pvalue = sapply(unique(pathways$Level2[pathways$Level1 == 'Metabolism']), test.metabolism.pathways)
metabolism.pvalue

test.energy.pathways = function(p){
  ctab = table(names(p.values) %in% 
                 pathways$KO[pathways$Level1 == 'Metabolism' & 
                               pathways$Level2 == 'Energy Metabolism' &
                               pathways$Level3 == p], 
               p.values<0.05)
  fisher.test(ctab)$p.value 
}
energy.pvalue = sapply(unique(pathways$Level3[pathways$Level1 == 'Metabolism' & pathways$Level2 == 'Energy Metabolism']), test.energy.pathways)
energy.pvalue

length(pathways$KO[pathways$Level3 == 'Photosynthesis'])

###
# Cecal
###
cv = function(x) abs(sd(x)/mean(x))
cecal = subset_samples(genes.deseq, Location == 'cecal')

cvs = apply(otu_table(cecal), 1, function(x) cv(x))
plot(ecdf(cvs))
plot(cvs, log='y')

cecal.var = subset_taxa(cecal, cvs>1)

Treatment = sample_data(cecal.var)$Treatment
x = otu_table(cecal.var)[1,]
summary(aov(c(x)~Treatment))


cecal.p = apply(otu_table(cecal.var), 1, function(x) summary(aov(c(x)~Treatment))[[1]][1,5])
cecal.q = p.adjust(cecal.p)
which(cecal.q < 0.05)
tax_table(cecal.var)["K00519",]

test.c.pathway.level1 = function(p){
  ctab = table(names(cecal.p) %in% pathways$KO[pathways$Level1 == p], cecal.p<0.05)
  fisher.test(ctab)$p.value
}
cecal.level1.pvalues = sapply(unique(pathways$Level1), test.c.pathway.level1)
cecal.level1.pvalues

## Metabolism is significant again
test.c.metabolism.pathways = function(p){
  ctab = table(names(cecal.p) %in% 
                 pathways$KO[pathways$Level1 == 'Metabolism' & pathways$Level2 == p], 
               cecal.p<0.05)
  fisher.test(ctab)$p.value
}
cecal.metabolism.pvalue = sapply(unique(pathways$Level2[pathways$Level1 == 'Metabolism']), test.c.metabolism.pathways)
cecal.metabolism.pvalue


