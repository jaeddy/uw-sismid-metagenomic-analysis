## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library("ape")
library("ggplot2")
library("ggtree")

## ------------------------------------------------------------------------
file <- system.file("extdata/BEAST", "beast_mcc.tree", package="ggtree")
beast <- read.beast(file)
beast

## ------------------------------------------------------------------------
get.fields(beast)

## ----warning=FALSE, fig.width=10, fig.height=10--------------------------
ggtree(beast, ndigits=2, branch.length = 'none') + geom_text(aes(x=branch, label=length_0.95_HPD), vjust=-.5, color='firebrick')

## ----warning=FALSE, fig.width=10, fig.height=10--------------------------
ggtree(beast) + geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=2)

## ----warning=FALSE, fig.width=10, fig.height=10--------------------------
ggtree(beast, branch.length = 'rate') + geom_range(range='rate_0.95_HPD', color='red', alpha=.6, size=2)

## ------------------------------------------------------------------------
beast_data <- fortify(beast)
head(beast_data)

## ----fig.width=12, fig.height=10, warning=FALSE, fig.align="center"------
brstfile <- system.file("extdata/PAML_Baseml", "rst", package="ggtree")
brst <- read.paml_rst(brstfile)
brst
p <- ggtree(brst) +
    geom_text(aes(x=branch, label=marginal_AA_subs), vjust=-.5, color='steelblue') +
    theme_tree2()
print(p)

## ------------------------------------------------------------------------
crstfile <- system.file("extdata/PAML_Codeml", "rst", package="ggtree")
crst <- read.paml_rst(crstfile)
crst

## ----fig.width=12, fig.height=10, warning=FALSE, fig.align="center"------
p %<% crst

## ------------------------------------------------------------------------
mlcfile <- system.file("extdata/PAML_Codeml", "mlc", package="ggtree")
mlc <- read.codeml_mlc(mlcfile)
mlc

## ----fig.width=8, fig.height=8, warning=FALSE, fig.align="center"--------
ggtree(mlc) + geom_text(aes(x=branch, label=dN_vs_dS), color='blue', vjust=-.2)

## ------------------------------------------------------------------------
get.fields(mlc)

## ----fig.width=8, width=60, warning=FALSE, fig.align="center"------------
ggtree(mlc, branch.length = "dN_vs_dS", aes(color=dN_vs_dS)) +
    scale_color_continuous(name='dN/dS', limits=c(0, 1.5),
                           oob=scales::squish, low="darkgreen", high="red")+
    theme_tree2(legend.position=c(.9, .5))

## ------------------------------------------------------------------------
ml <- read.codeml(crstfile, mlcfile)
ml
head(fortify(ml))

## ----warning=FALSE-------------------------------------------------------
nwk <- system.file("extdata/HYPHY", "labelledtree.tree", package="ggtree")
ancseq <- system.file("extdata/HYPHY", "ancseq.nex", package="ggtree")
tipfas <- system.file("extdata", "pa.fas", package="ggtree")
hy <- read.hyphy(nwk, ancseq, tipfas)
hy

## ----fig.width=16, fig.height=10, warning=FALSE, fig.align="center"------
ggtree(hy) + geom_text(aes(x=branch, label=AA_subs), vjust=-.5)

## ----fig.width=4, fig.height=6, width=60, warning=FALSE, fig.align="center"----
r8s <- read.r8s(system.file("extdata/r8s", "H3_r8s_output.log", package="ggtree"))
ggtree(r8s, branch.length="TREE", mrsd="2014-01-01") + theme_tree2()

## ----fig.width=16, fig.height=10, width=60, warning=FALSE, fig.align="center"----
ggtree(get.tree(r8s), aes(color=.id)) + facet_wrap(~.id, scales="free_x")

## ----fig.width=12, fig.height=10, width=60, warning=FALSE, fig.align="center"----
raxml_file <- system.file("extdata/RAxML", "RAxML_bipartitionsBranchLabels.H3", package="ggtree")
raxml <- read.raxml(raxml_file)
ggtree(raxml) + geom_label(aes(label=bootstrap, fill=bootstrap)) + geom_tiplab() +
    scale_fill_continuous(low='darkgreen', high='red') + theme_tree2(legend.position='right')

## ------------------------------------------------------------------------
nhxfile <- system.file("extdata", "ADH.nhx", package="ggtree")
nhx <- read.nhx(nhxfile)
ggtree(nhx) + geom_tiplab() + geom_point(aes(color=S), size=5, alpha=.5) + 
    theme(legend.position="right") +
    geom_text(aes(label=branch.length, x=branch), vjust=-.5) + 
    xlim(NA, 0.3)

## ------------------------------------------------------------------------
phyfile <- system.file("extdata", "sample.phy", package="ggtree")
phylip <- read.phylip(phyfile)
phylip
ggtree(phylip) + geom_tiplab()

## ----fig.width=10, fig.height=6------------------------------------------
msaplot(phylip, offset=1)

## ------------------------------------------------------------------------
jpf <- system.file("extdata/sample.jplace",  package="ggtree")
jp <- read.jplace(jpf)
print(jp)

## ------------------------------------------------------------------------
## get only best hit
get.placements(jp, by="best")
## get all placement
get.placements(jp, by="all")

## ------------------------------------------------------------------------
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)

rst_file <- system.file("examples/rst", package="ggtree")
mlc_file <- system.file("examples/mlc", package="ggtree")
codeml_tree <- read.codeml(rst_file, mlc_file)

merged_tree <- merge_tree(beast_tree, codeml_tree)

merged_tree
head(fortify(merged_tree))

## ----fig.width=20, fig.height=26, warning=FALSE--------------------------
ggtree(merged_tree, aes(color=dN_vs_dS), mrsd="2013-01-01", ndigits = 3) +
    geom_text(aes(label=posterior), vjust=.1, hjust=-.1, size=5, color="black") + 
    scale_color_continuous(name='dN/dS', limits=c(0, 1.5),
                           oob=scales::squish, low="green", high="red")+
    theme_tree2(legend.position="right")

## ------------------------------------------------------------------------
tree <- system.file("extdata", "pa.nwk", package="ggtree")
data <- read.csv(system.file("extdata", "pa_subs.csv", package="ggtree"), stringsAsFactor=FALSE)
print(tree)
head(data)

## ------------------------------------------------------------------------
outfile <- tempfile()
write.jplace(tree, data, outfile)

## ------------------------------------------------------------------------
jp <- read.jplace(outfile)
print(jp)

## ----fig.width=12, fig.height=12, warning=FALSE, fig.align="center"------
ggtree(jp) + 
    geom_text(aes(x=branch, label=subs), color="purple", vjust=-1, size=3) + 
    geom_text(aes(label=gc), color="steelblue", hjust=-.6, size=3) +
    geom_tiplab()

