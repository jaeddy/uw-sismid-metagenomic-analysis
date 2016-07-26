## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library("ape")
library("ggplot2")
library("ggtree")

## ------------------------------------------------------------------------
library("ggtree")
nwk <- system.file("extdata", "sample.nwk", package="ggtree")
tree <- read.tree(nwk)

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree, color="firebrick", size=1, linetype="dotted")

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree, ladderize=FALSE)

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree, branch.length="none")

## ----fig.height=4, fig.width=4, fig.align="center"-----------------------
ggtree(tree) + ggtitle("(Phylogram) rectangular layout")

## ----fig.height=4, fig.width=4, fig.align="center"-----------------------
ggtree(tree, layout="slanted") + ggtitle("(Phylogram) slanted layout")

## ----fig.height=5, fig.width=5, fig.align="center"-----------------------
ggtree(tree, layout="circular") + ggtitle("(Phylogram) circular layout")

## ----fig.height=5, fig.width=5, fig.align="center"-----------------------
ggtree(tree, layout="fan", open.angle=180) + ggtitle("(Phylogram) circular layout")

## ----fig.height=4, fig.width=4, fig.align="center"-----------------------
ggtree(tree, branch.length='none') + ggtitle("(Cladogram) rectangular layout")

## ----fig.height=4, fig.width=4, fig.align="center"-----------------------
ggtree(tree, layout="slanted", branch.length='none') + ggtitle("(Cladogram) slanted layout")

## ----fig.height=5, fig.width=5, fig.align="center"-----------------------
ggtree(tree, layout="circular", branch.length="none") + ggtitle("(Cladogram) circular layout")

## ----fig.height=5, fig.width=5, fig.align="center"-----------------------
ggtree(tree, layout="fan", open.angle=180, branch.length="none") + ggtitle("(Cladogram) circular layout")

## ----fig.height=4, fig.width=4, fig.align="center"-----------------------
ggtree(tree, layout="unrooted") + ggtitle("unrooted layout")

## ----fig.width=9, fig.height=9, fig.align="center"-----------------------
tree2d <- read.beast(system.file("extdata", "twoD.tree", package="ggtree"))
ggtree(tree2d, mrsd = "2014-05-01") + theme_tree2()

## ----fig.width=9, fig.height=4, fig.align="center"-----------------------
ggtree(tree2d, mrsd = "2014-05-01",
       yscale="NGS", yscale_mapping=c(N2=2, N3=3, N4=4, N5=5, N6=6, N7=7)) +
           theme_classic() + theme(axis.line.x=element_line(), axis.line.y=element_line()) +
               theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
                     panel.grid.major.y=element_blank()) +
                         scale_y_continuous(labels=paste0("N", 2:7))

## ----fig.width=4, fig.height=4, fig.align="center"-----------------------
ggtree(tree) + geom_treescale()

## ----fig.width=4, fig.height=4, fig.align="center"-----------------------
ggtree(tree)+geom_treescale(x=0, y=12, width=6, color='red')

## ----fig.width=4, fig.height=4, fig.align="center"-----------------------
ggtree(tree)+geom_treescale(fontsize=8, linesize=2, offset=-1)

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree) + theme_tree2()

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree)+geom_point(aes(shape=isTip, color=isTip), size=3)

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
p <- ggtree(tree) + geom_nodepoint(color="#b5e521", alpha=1/4, size=10)
p + geom_tippoint(color="#FDAC4F", shape=8, size=3)

## ----fig.width=3, fig.height=3, warning=FALSE, fig.align="center"--------
p + geom_tiplab(size=3, color="purple")

## ----fig.width=6, fig.height=6, warning=FALSE, fig.align="center"--------
ggtree(tree, layout="circular") + geom_tiplab(aes(angle=angle), color='blue')

## ----fig.width=6, fig.height=6, warning=FALSE, fig.align="center"--------
ggtree(tree, layout="circular") + geom_tiplab2(color='blue')

## ----fig.width=4, fig.height=3, warning=FALSE, fig.align="center"--------
p + geom_tiplab(aes(x=branch), size=3, color="purple", vjust=-0.3)

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
p %<% rtree(50)

## ----fig.width=6, fig.height=3, fig.align="center"-----------------------
multiplot(
    ggtree(rtree(30), color="red") + theme_tree("steelblue"),
    ggtree(rtree(20), color="white") + theme_tree("black"),
    ncol=2)

## ----fig.width=12, fig.height=4------------------------------------------
trees <- lapply(c(10, 20, 40), rtree)
class(trees) <- "multiPhylo"
ggtree(trees) + facet_wrap(~.id, scale="free") + geom_tiplab()

## ----fig.width=20, fig.height=20-----------------------------------------
btrees <- read.tree(system.file("extdata/RAxML", "RAxML_bootstrap.H3", package="ggtree"))
ggtree(btrees) + facet_wrap(~.id, ncol=10)

## ------------------------------------------------------------------------
p <- ggtree(btrees, layout="rectangular",   color="lightblue", alpha=.3)

best_tree <- read.tree(system.file("extdata/RAxML", "RAxML_bipartitionsBranchLabels.H3", package="ggtree"))
df <- fortify(best_tree, branch.length = 'none')
p+geom_tree(data=df, color='firebrick')

