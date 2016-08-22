## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library("ape")
library("ggplot2")
library("ggtree")

## ----fig.width=10, fig.height=5------------------------------------------
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)
beast_tree
p1 <- ggtree(beast_tree, mrsd='2013-01-01') + theme_tree2() +
    ggtitle("Divergence time")
p2 <- ggtree(beast_tree, branch.length = 'rate') + theme_tree2() +
    ggtitle("Substitution rate")
multiplot(p1, p2, ncol=2)

## ----fig.width=10, fig.height=5------------------------------------------
mlcfile <- system.file("extdata/PAML_Codeml", "mlc", package="ggtree")
mlc_tree <- read.codeml_mlc(mlcfile)
p1 <- ggtree(mlc_tree) + theme_tree2() +
    ggtitle("nucleotide substitutions per codon")
p2 <- ggtree(mlc_tree, branch.length='dN_vs_dS') + theme_tree2() +
    ggtitle("dN/dS tree")
multiplot(p1, p2, ncol=2)

## ------------------------------------------------------------------------
beast_tree2 <- rescale_tree(beast_tree, branch.length = 'rate')
ggtree(beast_tree2) + theme_tree2()

## ----fig.width=18, fig.height=10, fig.align="center"---------------------
library("ape")
data(chiroptera)
library("ggtree")
gzoom(chiroptera, grep("Plecotus", chiroptera$tip.label))

## ----fig.width=18, fig.height=10, warning=FALSE--------------------------
groupInfo <- split(chiroptera$tip.label, gsub("_\\w+", "", chiroptera$tip.label))
chiroptera <- groupOTU(chiroptera, groupInfo)
p <- ggtree(chiroptera, aes(color=group)) + geom_tiplab() + xlim(NA, 23)
gzoom(p, grep("Plecotus", chiroptera$tip.label), xmax_adjust=2)

## ----fig.width=5, fig.height=5-------------------------------------------
ggtree(beast_tree, aes(color=rate)) +
    scale_color_continuous(low='darkgreen', high='red') +
    theme(legend.position="right")

## ------------------------------------------------------------------------
set.seed(2015-12-21)
tree = rtree(30)
p <- ggtree(tree) + xlim(NA, 6)

p+geom_cladelabel(node=45, label="test label") +
    geom_cladelabel(node=34, label="another clade") 

## ------------------------------------------------------------------------
p+geom_cladelabel(node=45, label="test label", align=TRUE, offset=.5) +
    geom_cladelabel(node=34, label="another clade", align=TRUE, offset=.5)

## ------------------------------------------------------------------------
p+geom_cladelabel(node=45, label="test label", align=T, color='red') +
    geom_cladelabel(node=34, label="another clade", align=T, color='blue')

## ------------------------------------------------------------------------
p+geom_cladelabel(node=45, label="test label", align=T, angle=270, hjust='center', offset.text=.5) +
    geom_cladelabel(node=34, label="another clade", align=T, angle=45) 

## ------------------------------------------------------------------------
p+geom_cladelabel(node=45, label="test label", align=T, angle=270, hjust='center', offset.text=.5, barsize=1.5) +
    geom_cladelabel(node=34, label="another clade", align=T, angle=45, fontsize=8) 

## ------------------------------------------------------------------------
p+ geom_cladelabel(node=34, label="another clade", align=T, geom='label', fill='lightblue') 

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
nwk <- system.file("extdata", "sample.nwk", package="ggtree")
tree <- read.tree(nwk)
ggtree(tree) + geom_hilight(node=21, fill="steelblue", alpha=.6) +
    geom_hilight(node=17, fill="darkgreen", alpha=.6)

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
ggtree(tree, layout="circular") + geom_hilight(node=21, fill="steelblue", alpha=.6) +
    geom_hilight(node=23, fill="darkgreen", alpha=.6)

## ----fig.width=4, fig.height=5, fig.align='center', warning=FALSE--------
ggtree(tree) + 
  geom_balance(node=16, fill='steelblue', color='white', alpha=0.6, extend=1) +
  geom_balance(node=19, fill='darkgreen', color='white', alpha=0.6, extend=1)

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
ggtree(tree) + geom_tiplab() + geom_strip('E', 'G', barsize=2, color='red') + geom_strip('F', 'L', barsize=2, color='blue')

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
ggtree(tree) + geom_tiplab() + geom_taxalink('A', 'E') + geom_taxalink('F', 'K', color='red', arrow=grid::arrow(length = grid::unit(0.02, "npc")))

## ----results='hide', message=FALSE---------------------------------------
library(ape)
data(woodmouse)
d <- dist.dna(woodmouse)
tr <- nj(d)
bp <- boot.phylo(tr, woodmouse, function(xx) nj(dist.dna(xx)))

## ----fig.width=6, fig.height=6, warning=FALSE, fig.align="center"--------
tree <- apeBoot(tr, bp)
ggtree(tree) + geom_label(aes(label=bootstrap)) + geom_tiplab()

## ----results='hide', message=FALSE---------------------------------------
library(phangorn)
treefile <- system.file("extdata", "pa.nwk", package="ggtree")
tre <- read.tree(treefile)
tipseqfile <- system.file("extdata", "pa.fas", package="ggtree")
tipseq <- read.phyDat(tipseqfile,format="fasta")
fit <- pml(tre, tipseq, k=4)
fit <- optim.pml(fit, optNni=FALSE, optBf=T, optQ=T,
                 optInv=T, optGamma=T, optEdge=TRUE,
                 optRooted=FALSE, model = "GTR")

## ----fig.width=12, fig.height=10, width=60, warning=FALSE, fig.align="center"----
phangorn <- phyPML(fit, type="ml")
ggtree(phangorn) + geom_text(aes(x=branch, label=AA_subs, vjust=-.5))

## ------------------------------------------------------------------------
nwk <- system.file("extdata", "sample.nwk", package="ggtree")
tree <- read.tree(nwk)
p <- ggtree(tree)

dd <- data.frame(taxa  = LETTERS[1:13], 
                 place = c(rep("GZ", 5), rep("HK", 3), rep("CZ", 4), NA),
                 value = round(abs(rnorm(13, mean=70, sd=10)), digits=1))
## you don't need to order the data
## data was reshuffled just for demonstration
dd <- dd[sample(1:13, 13), ]
row.names(dd) <- NULL

## ----eval=FALSE----------------------------------------------------------
#  print(dd)

## ----echo=FALSE, results='asis'------------------------------------------
knitr::kable(dd)

## ----fig.width=6, fig.height=5, warning=FALSE, fig.align="center"--------
p <- p %<+% dd + geom_tiplab(aes(color=place)) + 
       geom_tippoint(aes(size=value, shape=place, color=place), alpha=0.25)
p+theme(legend.position="right")

## ----fig.width=6, fig.height=5, warning=FALSE, fig.align="center"--------
p + geom_text(aes(color=place, label=place), hjust=1, vjust=-0.4, size=3) +
    geom_text(aes(color=place, label=value), hjust=1, vjust=1.4, size=3)

## ----fig.width=6, fig.height=5, warning=FALSE, fig.align="center"--------
dd2 <- dd[, -1]
rownames(dd2) <- dd[,1]
require(phylobase)
tr2 <- phylo4d(tree, dd2)
ggtree(tr2) + geom_tiplab(aes(color=place)) + 
    geom_tippoint(aes(size=value, shape=place, color=place), alpha=0.25)

