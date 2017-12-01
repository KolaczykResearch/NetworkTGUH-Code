# Simulation of barbell network
library(far)
library(calibrate)
library(igraph)
library(matrixcalc)
require(graphics)
library(parallel)
library(rgl)

ptm <- proc.time()

n1 = 10  # Size of the first complete graph
n2 = 10 # Size of the second complete graph

C = barbell(n1, n2)  # Generate barbell
n = n1 + n2
pos = n*(n-1)/2  # Number of pairs of node
density = sum(C)/pos/2  # Density of the network
k = 100 # Numbers of network to simulate
retain = 2  # Number of the detail coefficients to retain

C.array = array(rep(as.matrix(C),k),c(n, n, k))  # Generate k copy of barbells
C.list = vector(mode = 'list', length = k)
C.list <- lapply(C.list, function(x) C.array[, , 1])  # Change the data type to list

q = 0.01  # Probability of committing type I error
p = (1-density)/density*q  # Probability of committing type II error

no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
result.list <- parLapply(cl, C.list, denoiseExample, p, q, retain, 0)
stopCluster(cl)

result = lapply(C.list, tguhBarbell, p, q)

detail = lapply(result, unlistCompression)
temp = matrix(unlist(detail), nrow=100, byrow=TRUE)
mean.details = apply(temp,2,mean)

C = barbell(n1, n2)  # Generate barbell
trueB = uh.bu.net.nonrem(C, column.norm, 2)
test = uh.bu.net.nonrem.mult(C, column.norm, 2, 0.01)


truebarbell = unlist(lapply(trueB$detail.list, function(y) column.norm(as.matrix(y),p=2)))

plotCompression(truebarbell, mean.details)

proc.time() - ptm

# Denoise of the DTI network
library(rgl)
library(igraph)
library(Matrix)

setwd("~/Documents/BU/Research/Fast bottom-up Decompostion Algorithm on Network/BrainNetwork Data")
C = read.csv('Connectivity.csv', header=FALSE)
S = read.csv('PowerSpectra.csv', header=FALSE)
xyz = read.csv('Coordinates.csv', header=FALSE)

C = sym.by.avg(as.matrix(C))$mat # make the adjacency matrix a symmetric one

# Make C a symmetric matrix
C = sym.by.avg(as.matrix(C))$mat # make the adjacency matrix a symmetric one

# Pre-process to remove some noise in connectivities
logvalue = -10
sum(as.vector(log(C))>logvalue) # check number of edges

g = graph_from_adjacency_matrix(as.matrix(log(C)>logvalue)*C, mode="undirected", weighted=TRUE) # make the data.frame a connectivity network
is.connected(g)  # test if the graph is connected

g = normalize.signal(g, S[,2]) # attach normalized signal over the network (S[,1]=theta, S[,2]=alpha, S[,3]=beta, S[,4]=gamma)

CC = decompose(g, max.comps=1) # extract the leading giant connected component if g is not a connected graph

# Extract the (simple) adjacency matrix of g 
# C = as.matrix(as_adjacency_matrix(CC[[1]], attr="weight"))  # weighted version
C = as.matrix(as_adjacency_matrix(CC[[1]]))                   # unweighted version
diag(C) = V(CC[[1]])$signal

# trans.wth.ns = uh.bu.net.nonrem.mult(C, column.norm, 2, 0.01)      # tail greedy version
trans.wth.ns = uh.bu.net.nonrem(C, column.norm, 2)                   # greedy version

denoised = uh.bu.net.inv.th(trans.wth.ns, NA) # th = NA means the theoretical threshold will be used
sig.den = diag(denoised[[1]])

series = c(V(g)$signal, sig.den)
Col = computeColor(series, 1000)
Col[,1] = as.numeric(levels(Col[,1]))[Col[,1]]

rglplot(g, layout=as.matrix(xyz), vertex.label=NA, edge.label=NA, vertex.size=V(g)$signal, vertex.color=Col[match(round(V(g)$signal, 3), Col[, 1]),2])
title3d(xlab="x", ylab="y", zlab="z")
rgl.bbox(xlen=max(xyz), ylen=max(xyz), zlen=max(xyz))

rglplot(g, layout=as.matrix(xyz), vertex.label=NA, edge.label=NA, vertex.size=sig.den, vertex.color=Col[match(round(sig.den, 3), Col[, 1]),2])
title3d(xlab="x", ylab="y", zlab="z")
rgl.bbox(xlen=max(xyz), ylen=max(xyz), zlen=max(xyz))
