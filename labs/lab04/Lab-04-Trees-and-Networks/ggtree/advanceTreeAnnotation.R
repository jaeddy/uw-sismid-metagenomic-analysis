## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library("ape")
library("ggplot2")
library("ggtree")

## ----fig.width=20, fig.height=16, fig.align="center"---------------------
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)

genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
p <- ggtree(beast_tree, mrsd="2013-01-01") + geom_treescale(x=2008, y=1)
p <- p + geom_tiplab(size=3)
gheatmap(p, genotype, offset = 2, width=0.5)

## ----fig.width=20, fig.height=16, fig.align="center", warning=FALSE------
p <- ggtree(beast_tree, mrsd="2013-01-01") + geom_tiplab(size=3, align=TRUE) + theme_tree2()
pp <- (p + scale_y_continuous(expand=c(0, 0.3))) %>%
    gheatmap(genotype, offset=4, width=0.5, colnames=FALSE) %>%
        scale_x_ggtree()
pp + theme(legend.position="right")

## ----fig.width=8, fig.height=12, fig.align='center'----------------------
fasta <- system.file("examples/FluA_H3_AA.fas", package="ggtree")
msaplot(ggtree(beast_tree), fasta) 

## ----fig.width=16, fig.height=16, fig.align='center'---------------------
msaplot(ggtree(beast_tree), fasta, window=c(150, 200)) + coord_polar(theta='y')

## ------------------------------------------------------------------------
set.seed(2015-12-31)
tr <- rtree(15)
p <- ggtree(tr)

a <- runif(14, 0, 0.33)
b <- runif(14, 0, 0.33)
c <- runif(14, 0, 0.33)
d <- 1 - a - b - c
dat <- data.frame(a=a, b=b, c=c, d=d)
## input data should have a column of `node` that store the node number
dat$node <- 15+1:14

## cols parameter indicate which columns store stats (a, b, c and d in this example)
bars <- nodebar(dat, cols=1:4)

inset(p, bars)

## ------------------------------------------------------------------------
inset(p, bars, width=.03, height=.06)

## ------------------------------------------------------------------------
bars2 <- nodebar(dat, cols=1:4, position='dodge',
                 color=c(a='blue', b='red', c='green', d='cyan'))
p2 <- inset(p, bars2, x='branch', width=.03, vjust=-.3)
print(p2)

## ------------------------------------------------------------------------
pies <- nodepie(dat, cols=1:4, alpha=.6)
inset(p, pies)

## ------------------------------------------------------------------------
inset(p, pies, hjust=-.06)

## ------------------------------------------------------------------------
pies_and_bars <- bars2
pies_and_bars[9:14] <- pies[9:14]
inset(p, pies_and_bars)

## ------------------------------------------------------------------------
d <- lapply(1:15, rnorm, n=100)
ylim <- range(unlist(d))
bx <- lapply(d, function(y) {
    dd <- data.frame(y=y)
    ggplot(dd, aes(x=1, y=y))+geom_boxplot() + ylim(ylim) + theme_inset()
})
names(bx) <- 1:15
inset(p, bx, width=.03, height=.1, hjust=-.05)

## ----fig.width=10, fig.height=7------------------------------------------
p2 <- inset(p, bars2, x='branch', width=.03, vjust=-.4)
p2 <- inset(p2, pies, x='branch', vjust=.4)
bx2 <- lapply(bx, function(g) g+coord_flip())
inset(p2, bx2, width=.2, height=.03, vjust=.04, hjust=p2$data$x[1:15]-4) + xlim(NA, 4.5)

## ----eval=FALSE----------------------------------------------------------
#  imgfile <- tempfile(, fileext=".png")
#  download.file("https://avatars1.githubusercontent.com/u/626539?v=3&u=e731426406dd3f45a73d96dd604bc45ae2e7c36f&s=140", destfile=imgfile, mode='wb')
#  img <- list(imgfile, imgfile)
#  names(img) <- c("18", "22")
#  inset(p, img)

## ----warning=F, fig.width=10, fig.height=6-------------------------------
tr <- rtree(30)
df <- fortify(tr)
df$tipstats <- NA
d1 <- df
d2 <- df
d2$tipstats[d2$isTip] <- abs(rnorm(30))
d1$panel <- 'Tree'
d2$panel <- 'Stats'
d1$panel <- factor(d1$panel, levels=c("Tree", "Stats"))
d2$panel <- factor(d2$panel, levels=c("Tree", "Stats"))

p <- ggplot(mapping=aes(x=x, y=y)) + facet_grid(.~panel, scale="free_x") + theme_tree2()
p+geom_tree(data=d1) + geom_point(data=d2, aes(x=tipstats)) 

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE, eval=FALSE----
#  pp <- ggtree(tree) %>% phylopic("79ad5f09-cf21-4c89-8e7d-0c82a00ce728", color="steelblue", alpha = .3)
#  print(pp)

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE, eval=FALSE----
#  pp %>% phylopic("67382184-5135-4faa-8e98-eadff02c3e8a", color="#86B875", alpha=.8, node=4) %>%
#       phylopic("d3563b54-780f-4711-a49a-7ea051e9dacc", color="darkcyan", alpha=.8, node=17, width=.2)

