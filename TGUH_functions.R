barbell <- function(n1, n2){
# simulate barbell network of size n1 + n2
  m1 = matrix(1, n1, n1)
  m2 = matrix(1, n2, n2)
  
  block.mat = bdiag(m1, m2)
  block.mat[n1, n1+1] = 1
  block.mat[n1+1, n1] = 1
  
  return(block.mat)
}

diagSparse <- function(n1, n2){
# Build a block diagonal matrix given several building block matrices.
  m1 = matrix(1, n1, n1)
  m2 = matrix(1, n1, n2)
  
  return(bdiag(m1, m2))
}

one <- function(n, p = n){
# generate matrix of ones
  return(matrix(rep(1,n*n),n,n))
}

tguhBarbell <- function(mat, p, q){
  # tguh transform of noisy barbell network
  row = nrow(mat)
  col = ncol(mat)
  count = 0
  
  L2.mat = NULL
  L2.W.hat = NULL
  L2.noise = NULL
  
  typeI = apply(mat,1,function(x) runif(row))
  typeII = apply(mat,1,function(x) runif(row))
  
  U = round(upper.tri(mat))
  
  P = ((mat*typeII)>p)*1 #Gen Type II err: For those =1, entry with zero stays zero, entry with one become become zero with prob p
  Q = ((mat+typeI)>(1-q))*1 #Gen Type I err: For those =0, any entry > 1-p becomes 1, otherwise stay 0
  
  W.kn = matrix(rep(1,row*row),row,row) - diag(row)
  upper = as.matrix((P+Q-mat)*U)
  
  W.hat = t(upper)+upper
  diag(W.hat) = 1
  
  if(!all(W.hat == mat)){
    count = count + 1
  }
  trs.mat = uh.bu.net.nonrem(mat, column.norm, 2)
  trs.W.hat = uh.bu.net.nonrem(W.hat, column.norm, 2)
  
  return(list(trs = trs.mat, trs.W = trs.W.hat, W.hat=W.hat, count))
}

unlistCompression <- function(x){ # original fun2
  # unlist the compression of simulations
  unlist(lapply(x$trs.W$detail.list, function(y) column.norm(as.matrix(y),p=2)))
}

computeColor <- function(series, n){
  # compute color for the rgl plot
  
  value.max = max(series)
  value.min = min(series)
  n = round(value.max-value.min, 3)*1000
  rbPal <- colorRampPalette(c('red','blue'))
  
  
  sequence = seq(value.min, value.max, (value.max-value.min)/(n-1))
  Col <- rbPal(n)[as.numeric(cut(sequence,breaks=n))]
  df = as.data.frame(cbind(round(sequence, 3),Col))
  
  return(df)
}

denoiseExample <- function(mat, p, q, retain, count){
  
  vec.norm <- function(x, p=2) {
    
    sum(abs(x)^p)^(1/p)
  }
  
  column.norm <- function(A, p=2) {
    #	Vector norm (eg L1, p=1, or L2, p=2) applied to matrices, for measuring
    #	similarity between the columns of the adjacency matrix. Use as column.norm
    # 	in uh.bu.net.
    
    apply(A, p, vec.norm, p)
  }
  
  uh.bu.net.inv.sm <- function(a, retain) {
    #	Inverse transform with "smoothing": all but the 'retain' largest detail columns get set
    #	to zero before the inverse transform is taken. This is what I used to produce the
    # denoised matrices for the slides.
    
    n <- dim(a$A)[1]
    for (s in (n-1):1) {
      
      eating.up <- min(a$decomp.hist[s, 1:2])
      eaten.up <- max(a$decomp.hist[s, 1:2])
      
      h1 <- a$decomp.hist[s, 3]
      
      mat <- matrix(c(sqrt(1-h1^2), h1, h1, -sqrt(1-h1^2)), 2, 2)
      
      if (s > (n-1-retain)) {
        a$A[c(eating.up, eaten.up),] <- mat %*% a$A[c(eating.up, eaten.up),]
        a$A[,c(eating.up, eaten.up)] <- a$A[,c(eating.up, eaten.up)] %*% mat
      }
      else {
        a$A[eaten.up,] <- 0
        a$A[,eaten.up] <- 0
        a$A[c(eating.up, eaten.up),] <- mat %*% a$A[c(eating.up, eaten.up),]
        a$A[,c(eating.up, eaten.up)] <- a$A[,c(eating.up, eaten.up)] %*% mat
      }
    }
    a$A
  }
  
  uh.bu.net.nonrem <- function(A, column.norm, p) {
    #	  This is the same as uh.bu.net (=the two algorithms lead to identical results), the only difference being
    #   in the implementation: in this one, the matrix A does not decrease in size as the algorithm progresses,
    #   which is useful for the inverse transform (below).
    
    n <- dim(A)[1]
    B <- A
    B[lower.tri(B, diag=T)] <- 0
    noe <- sum(B == 1)
    weights <- rep(1, n)
    edges <-  matrix(which(B == 1, arr.ind=T), noe, 2)
    decomp.hist <- matrix(0, n-1, 3)
    cluster.hist <- t(matrix(1:n, n, n))
    active.rows <- rep(1, n)
    detail.list <- list()
    
    for (s in 1:(n-1)) {
      #	print(s)
      a <- weights[edges[,1]]
      b <- weights[edges[,2]]
      
      h1 <- 1/sqrt(1 + (a/b)^2)
      h2 <- -1/sqrt(1 + (b/a)^2)
      
      l1 <- -h2
      l2 <- h1
      
      details <- t(t(A[,edges[,1]]) * h1 + t(A[,edges[,2]]) * h2)
      
      details.min.ind <- which.min(column.norm(details * active.rows, p))
      
      if (length(details.min.ind) > 1) details.min.ind <- sample(details.min.ind, 1)
      
      smooth.at.min <- A[,edges[details.min.ind,1]] * l1[details.min.ind] + A[,edges[details.min.ind,2]] * l2[details.min.ind]
      detail.at.min <- A[,edges[details.min.ind,1]] * h1[details.min.ind] + A[,edges[details.min.ind,2]] * h2[details.min.ind]
      
      sm.weight.at.min <- l1[details.min.ind] * weights[edges[details.min.ind,1]] + 
        l2[details.min.ind] * weights[edges[details.min.ind,2]]
      
      decomp.hist[s, 1:2] <- edges[details.min.ind,]
      decomp.hist[s,3] <- h1[details.min.ind]
      
      detail.list[[s]] <- details[,details.min.ind]		
      
      eating.up <- min(edges[details.min.ind,])
      eaten.up <- max(edges[details.min.ind,])
      
      A[,eating.up] <- smooth.at.min
      A[,eaten.up] <- detail.at.min
      smooth.at.min <- A[eating.up,] * l1[details.min.ind] + A[eaten.up,] * l2[details.min.ind]
      detail.at.min <- A[eating.up,] * h1[details.min.ind] + A[eaten.up,] * h2[details.min.ind]
      
      A[eating.up,] <- smooth.at.min
      A[eaten.up,] <- detail.at.min
      
      active.rows[eaten.up] <- 0
      
      weights[eating.up] <- sm.weight.at.min
      
      ind <- which(edges == eaten.up)
      edges[ind] <- eating.up
      d <- dim(edges)[1]
      rows.changed <- union(ind[ind > d] - d, ind[ind <= d])
      rows.changed.equal <- rows.changed[edges[rows.changed,1] == edges[rows.changed,2]]		
      edges <- edges[-rows.changed.equal,,drop=F]
      edges <- t(apply(edges, 1, sort))
      
      #		print(eaten.up)
      #		print(which(cluster.hist[s,] == eaten.up))
      
      cluster.hist[s, which(cluster.hist[s,] == eaten.up)] <- eating.up
      cluster.hist[s+1,] <- cluster.hist[s,]
      
      if (dim(edges)[1] == 0) 
        break
      if (length(edges) == 2) 
        edges <- matrix(edges, 1, 2)
      
    }
    list(A=A, decomp.hist=decomp.hist, weights=weights, detail.list=detail.list, cluster.hist=cluster.hist)
  }
  #p: probability of removing a existing edge (type II err)
  #q: probability of adding a non-existing edge (type I err)
  
  row = nrow(mat)
  col = ncol(mat)
  
  L2.mat = NULL
  L2.W.hat = NULL
  L2.noise = NULL
  
  typeI = apply(mat,1,function(x) runif(row))
  typeII = apply(mat,1,function(x) runif(row))
  
  U = round(upper.tri(mat))
  
  P = ((mat*typeII)>p)*1 #Gen Type II err: For those =1, entry with zero stays zero, entry with one become become zero with prob p
  Q = ((mat+typeI)>(1-q))*1 #Gen Type I err: For those =0, any entry > 1-p becomes 1, otherwise stay 0
  
  W.kn = matrix(rep(1,row*row),row,row) - diag(row)
  upper = as.matrix((P+Q-mat)*U)
  
  W.hat = t(upper)+upper
  diag(W.hat) = 1
  
  if(!all(W.hat == mat)){
    count = count + 1
  }
  trs.mat = uh.bu.net.nonrem(mat, column.norm, 2)
  trs.W.hat = uh.bu.net.nonrem(W.hat, column.norm, 2)
  
  
  
  for(i in retain:(row)){
    de.mat = uh.bu.net.inv.sm(trs.mat, i)
    de.W.hat = uh.bu.net.inv.sm(trs.W.hat, i)
    
    tol = 1e-5
    de.mat = as.matrix((abs(de.mat) > tol) * de.mat)
    de.W.hat = as.matrix((abs(de.W.hat) > tol) * de.W.hat)
    
    L2.mat[i-1] = norm(de.mat - W.hat)
    L2.W.hat[i-1] = norm(de.W.hat - W.hat)
  }
  
  L2.noise = norm(W.hat - mat, type='f')
  
  list(L2.mat, L2.W.hat, rep(L2.noise,row-retain+1), rep(count,row-retain+1), trs.W.hat$decomp.hist)
  
  #  return(t(upper)+upper) # W.true + add(P - W.true) - remove(W.true - Q); #W.kn remove all self loops
}

vec.norm <- function(x, p) {
  
  sum(abs(x)^p)^(1/p)
}

row.norm <- function(A, p) {
  #	Vector norm (eg L1, p=1, or L2, p=2) applied to matrices, for measuring
  #	similarity between the columns of the adjacency matrix. 
  
  apply(A, 1, vec.norm, p)
}

column.norm <- function(A, p) {
  #	Vector norm (eg L1, p=1, or L2, p=2) applied to matrices, for measuring
  #	similarity between the columns of the adjacency matrix. Use as column.norm
  # 	in uh.bu.net.
  
  apply(A, 2, vec.norm, p)
}

uh.bu.net.inv.sm <- function(a, retain) {
  #	Inverse transform with "smoothing": all but the 'retain' largest detail columns get set
  #	to zero before the inverse transform is taken.
  
  n <- dim(a$A)[1]
  for (s in (n-1):1) {
    
    eating.up <- min(a$decomp.hist[s, 1:2])
    eaten.up <- max(a$decomp.hist[s, 1:2])
    
    h1 <- a$decomp.hist[s, 3]
    
    mat <- matrix(c(sqrt(1-h1^2), h1, h1, -sqrt(1-h1^2)), 2, 2)
    
    if (s > (n-1-retain)) {
      a$A[c(eating.up, eaten.up),] <- mat %*% a$A[c(eating.up, eaten.up),]
      a$A[,c(eating.up, eaten.up)] <- a$A[,c(eating.up, eaten.up)] %*% mat
    }
    else {
      a$A[eaten.up,] <- 0
      a$A[,eaten.up] <- 0
      a$A[c(eating.up, eaten.up),] <- mat %*% a$A[c(eating.up, eaten.up),]
      a$A[,c(eating.up, eaten.up)] <- a$A[,c(eating.up, eaten.up)] %*% mat
    }
  }
  a$A
}

uh.bu.net.inv.th <- function(a, th) {
  #	Similar to function uh.bu.net.inv.sm.
  # th is the threshold provided for the inverse transform.
  # If th is NA, then the theorectical value will be used.
  
  n <- dim(a$A)[1]
  if(is.na(th))
    th = 2*log(n)
  count = 1
  
  for (s in (n-1):1) {
    
    eating.up <- min(a$decomp.hist[s, 1:2])
    eaten.up <- max(a$decomp.hist[s, 1:2])
    
    h1 <- a$decomp.hist[s, 3]
    
    mat <- matrix(c(sqrt(1-h1^2), h1, h1, -sqrt(1-h1^2)), 2, 2)
    
    if (all(row.norm(a$A[c(eating.up, eaten.up),],2)>th)) {
      a$A[c(eating.up, eaten.up),] <- mat %*% a$A[c(eating.up, eaten.up),]
      a$A[,c(eating.up, eaten.up)] <- a$A[,c(eating.up, eaten.up)] %*% mat
      count = count+1
    }
    else {
      a$A[eaten.up,] <- 0
      a$A[,eaten.up] <- 0
      a$A[c(eating.up, eaten.up),] <- mat %*% a$A[c(eating.up, eaten.up),]
      a$A[,c(eating.up, eaten.up)] <- a$A[,c(eating.up, eaten.up)] %*% mat
    }
  }
  list(a$A, count, th)
}

uh.bu.net.nonrem <- function(A, column.norm, p) {
  #   This is the greedy version of the TGUH where only one node if merged at each run
  
  n <- dim(A)[1]
  B <- A
  B[lower.tri(B, diag=T)] <- 0
  noe <- sum(B == 1)
  weights <- rep(1, n)
  edges <-  matrix(which(B == 1, arr.ind=T), noe, 2)
  decomp.hist <- matrix(0, n-1, 3)
  cluster.hist <- t(matrix(1:n, n, n))
  active.rows <- rep(1, n)
  detail.list <- list()
  
  for (s in 1:(n-1)) {
    #	print(s)
    a <- weights[edges[,1]]
    b <- weights[edges[,2]]
    
    h1 <- 1/sqrt(1 + (a/b)^2)
    h2 <- -1/sqrt(1 + (b/a)^2)
    
    l1 <- -h2
    l2 <- h1
    
    details <- t(t(A[,edges[,1]]) * h1 + t(A[,edges[,2]]) * h2)
    
    details.min.ind <- which.min(column.norm(details * active.rows, p))
    
    if (length(details.min.ind) > 1) details.min.ind <- sample(details.min.ind, 1)
    
    smooth.at.min <- A[,edges[details.min.ind,1]] * l1[details.min.ind] + A[,edges[details.min.ind,2]] * l2[details.min.ind]
    detail.at.min <- A[,edges[details.min.ind,1]] * h1[details.min.ind] + A[,edges[details.min.ind,2]] * h2[details.min.ind]
    
    sm.weight.at.min <- l1[details.min.ind] * weights[edges[details.min.ind,1]] + 
      l2[details.min.ind] * weights[edges[details.min.ind,2]]
    
    decomp.hist[s, 1:2] <- edges[details.min.ind,]
    decomp.hist[s,3] <- h1[details.min.ind]
    
    detail.list[[s]] <- details[,details.min.ind]		
    
    eating.up <- min(edges[details.min.ind,])
    eaten.up <- max(edges[details.min.ind,])
    
    A[,eating.up] <- smooth.at.min
    A[,eaten.up] <- detail.at.min
    smooth.at.min <- A[eating.up,] * l1[details.min.ind] + A[eaten.up,] * l2[details.min.ind]
    detail.at.min <- A[eating.up,] * h1[details.min.ind] + A[eaten.up,] * h2[details.min.ind]
    
    A[eating.up,] <- smooth.at.min
    A[eaten.up,] <- detail.at.min
    
    active.rows[eaten.up] <- 0
    
    weights[eating.up] <- sm.weight.at.min
    
    ind <- which(edges == eaten.up)
    edges[ind] <- eating.up
    d <- dim(edges)[1]
    rows.changed <- union(ind[ind > d] - d, ind[ind <= d])
    rows.changed.equal <- rows.changed[edges[rows.changed,1] == edges[rows.changed,2]]		
    edges <- edges[-rows.changed.equal,,drop=F]
    edges <- t(apply(edges, 1, sort))
    
    #		print(eaten.up)
    #		print(which(cluster.hist[s,] == eaten.up))
    
    cluster.hist[s, which(cluster.hist[s,] == eaten.up)] <- eating.up
    cluster.hist[s+1,] <- cluster.hist[s,]
    
    if (dim(edges)[1] == 0) 
      break
    if (length(edges) == 2) 
      edges <- matrix(edges, 1, 2)
    
  }
  list(A=A, decomp.hist=decomp.hist, weights=weights, detail.list=detail.list, cluster.hist=cluster.hist)
}

uh.bu.net.nonrem.mult <- function(A, column.norm, p, rho) {
  # This is the tguh version of the algorithm, where multiple node is allowed to merge at each step
  # Parameter rho determines the speedy/greedyness of he algorithm
  
  n <- dim(A)[1]
  B <- A
  B[lower.tri(B, diag=T)] <- 0
  #	noe <- sum(B == 1) # for simple graph
  noe <- sum(B != 0) # for weighted graph
  weights <- rep(1, n)
  edges <-  matrix(which(B != 0 , arr.ind=T), noe, 2)
  decomp.hist <- matrix(0, n-1, 3)
  cluster.hist <- t(matrix(1:n, n, n))
  active.rows <- rep(1, n)
  detail.list <- list()
  active.level <- rep(1, n)
  remain = n
  
  ind2 = 1
  level = 1
  
  while (dim(edges)[2] != 0) {
    max.size = ceiling(remain*rho)
    for(j in 1:max.size){
      #	print(s)
      a <- weights[edges[,1]]
      b <- weights[edges[,2]]
      
      h1 <- 1/sqrt(1 + (a/b)^2)
      h2 <- -1/sqrt(1 + (b/a)^2)
      
      l1 <- -h2
      l2 <- h1
      
      details <- t(t(A[,edges[,1]]) * h1 + t(A[,edges[,2]]) * h2) # compute details column by column
      
      details.min.ind <- which.min(column.norm(details* (active.rows *active.level), p)) # select the column with the lowest norm.
      
      if (length(details.min.ind) > 1) details.min.ind <- sample(details.min.ind, 1)  # random choose one column if multiple equal to each other
      
      # compute the meta node and the detail difference
      smooth.at.min <- A[,edges[details.min.ind,1]] * l1[details.min.ind] + A[,edges[details.min.ind,2]] * l2[details.min.ind]
      detail.at.min <- A[,edges[details.min.ind,1]] * h1[details.min.ind] + A[,edges[details.min.ind,2]] * h2[details.min.ind]
      
      sm.weight.at.min <- l1[details.min.ind] * weights[edges[details.min.ind,1]] + 
        l2[details.min.ind] * weights[edges[details.min.ind,2]]
      
      decomp.hist[ind2, 1:2] <- edges[details.min.ind,]
      decomp.hist[ind2, 3] <- h1[details.min.ind]
      
      detail.list[[ind2]] <- details[,details.min.ind]		
      
      eating.up <- min(edges[details.min.ind,])
      eaten.up <- max(edges[details.min.ind,])
      
      A[,eating.up] <- smooth.at.min
      A[,eaten.up] <- detail.at.min
      smooth.at.min <- A[eating.up,] * l1[details.min.ind] + A[eaten.up,] * l2[details.min.ind]
      detail.at.min <- A[eating.up,] * h1[details.min.ind] + A[eaten.up,] * h2[details.min.ind]
      
      A[eating.up,] <- smooth.at.min
      A[eaten.up,] <- detail.at.min
      
      active.rows[eaten.up] <- 0
      active.level[eating.up] <- 0
      
      weights[eating.up] <- sm.weight.at.min
      
      ind <- which(edges == eaten.up)
      edges[ind] <- eating.up
      d <- dim(edges)[1]
      rows.changed <- union(ind[ind > d] - d, ind[ind <= d])
      rows.changed.equal <- rows.changed[edges[rows.changed,1] == edges[rows.changed,2]]		
      edges <- edges[-rows.changed.equal,,drop=F]
      edges <- t(apply(edges, 1, sort))
      
      #		print(eaten.up)
      #		print(which(cluster.hist[s,] == eaten.up))
      
      cluster.hist[ind2, which(cluster.hist[ind2,] == eaten.up)] <- eating.up
      cluster.hist[ind2,] <- cluster.hist[ind2,]
      
      if (length(edges) == 2) 
        edges <- matrix(edges, 1, 2)
      ind2 = ind2+1
      remain = remain-1
    }
    active.level <- rep(1, n)
    level = level+1
    dim(edges)[2]
  }
  list(A=A, decomp.hist=decomp.hist, weights=weights, detail.list=detail.list, cluster.hist=cluster.hist, level=level)
}

plotCompression <- function(true, noise){
  
  true = true*(true>1)
  plot(rev(noise), type='l', ylab=NA, xlab=NA)
  
  title(ylab = "Detail coefficients (Avg)", xlab='Coefficients index', line=2, cex.lab=1)
  
  lines(rev(true), lty=3)
  legend("topright",legend=c("noisy", "noise free"), lty=c(1,3))
}

sym.by.avg <- function(A) {
  # make a asymmetric matrix a symmetric one
  # by taking the average of the two directions
  l = lower.tri(A)*A
  u = upper.tri(A)*A
  
  newupper = (u + t(l))/2
  diag(newupper) = diag(A)
  mat = forceSymmetric(newupper)
  
  list(mat = mat)
}

vec.inverse = function(i, A){
  # vectorized inverse transform
  
  uh.bu.net.inv.sm(A, i)
}

normalize.signal <- function(g, signal) {
  # normalize f function to unit standard deviations
  noe = length(E(g))
  
  V(g)$signal = signal
  sigma.hat = est.noise(g)
  V(g)$signal = V(g)$signal/sigma.hat
  
  #  h = V(g)$signal[head_of(g, 1:noe)]
  #  t = V(g)$signal[tail_of(g, 1:noe)]
  
  #  E(g)$weight = (h+t)/2
  return(g)
}

est.noise <- function(g) {
  
  # Estimates the standard deviation of iid Gaussian noise
  # in network Median Absolute Deviation 
  # (Davies and Kovac 2001, Local Extremes, Runs, Strings, and Multiresolution).
  noe = length(E(g))
  
  h = V(g)$signal[head_of(g, 1:noe)]
  t = V(g)$signal[tail_of(g, 1:noe)]
  
  d = h-t
  mad(d/sqrt(2))
}
