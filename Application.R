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


C = barbell(4,6)
g = graph_from_adjacency_matrix(C, mode="undirected", weighted=TRUE) # make the 
is.connected(g)  # test if the graph is connected
set.seed(1)
spn = c(rep(0,3),3,7,rep(10,5)) + rnorm(10)  # generate signal and noise
# Extract the (simple) adjacency matrix of g 
C = as.matrix(as_adjacency_matrix(g))

trans = uh.bu.net.nonrem.mult(C, column.norm, 2, 0.01)      # tail greedy version
denoised.signal = denoise.th(trans, spn, NA)                # th = NA means the theoretical threshold will be used
