library(devtools)
install_github('zdk123/SpiecEasi')
library(SpiecEasi)

data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))
d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))
X <- synth_comm_from_counts(amgut1.filt.cs, 
                            mar=2, 
                            distr='zinegbin', 
                            Sigma=Cor, n=n)

se.est <- spiec.easi(X, method='mb', 
                     lambda.min.ratio=1e-2, 
                     nlambda=15)

ig <- graph.adjacency(graph, mode='undirected')
# save layout for side-by-side plotting
g.coord <- layout.fruchterman.reingold(ig)
plot(ig, layout=g.coord, vertex.size=6, vertex.label=NA)
plot(graph.adjacency(se.est$refit, mode='undirected'),
     layout=g.coord, vertex.size=6, vertex.label=NA)

huge::huge.roc(se.est$path, graph)

SpiecEasi:::stars.pr(se.est$merge[[se.est$opt.index]], graph, ll=15)



### American Gut
se.mb.amgut <- 
  spiec.easi(amgut1.filt, method='mb', 
             lambda.min.ratio=1e-2, nlambda=20, 
             icov.select.params=list(rep.num=50))

se.gl.amgut <- 
  spiec.easi(amgut1.filt, method='glasso', 
             lambda.min.ratio=1e-2, nlambda=20, 
             icov.select.params=list(rep.num=50))

sparcc.amgut <- sparcc(amgut1.filt)
## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3

## Create igraph objects
ig.mb <- graph.adjacency(se.mb.amgut$refit, mode='undirected')
ig.gl <- graph.adjacency(se.gl.amgut$refit, mode='undirected')
ig.sparcc <- graph.adjacency(sparcc.graph, mode='undirected', diag=FALSE)

## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(amgut1.filt, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)
par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize,
      vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, 
     vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, 
     vertex.label=NA, main="sparcc")

dd.gl <- degree.distribution(ig.gl)
dd.mb <- degree.distribution(ig.mb)
dd.sparcc <- degree.distribution(ig.sparcc)

par(mfrow=c(1,1))
plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,.35), pch=19, type='b',
     ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b', pch=19)
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b', pch=19)
legend("topright", c("MB", "glasso", "sparcc"),
       col=c("forestgreen", "red", "black"), pch=19, lty=1)


