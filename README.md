
README
Xinyu Kang
xkang@bu.edu
December 2017

This is the TGUH algorithm that implemented the tail greedy unbalanced haar transformation of networks.

Multiscale network analysis through tail-greedy bottom-up approximation, with applications in neuroscience (2017)
Xinyu Kang, Piotr Fryzlewicz, Catherine Chu, Mark Kramer, Eric Kolaczyk
Asilomar Conference on Signals, Systems, and Computers (2017)



Applications 
=============================================================================

- Part 1: 	Simulation study that generates Figure 1 in the paper

- Part 2: 	An denoise example using the network TGUH.
		We denoise the signal over a barbell of size 4+6.

Functions 
=============================================================================

TGUH_functions.r : contains  functions used by the algorithm.

barbell			: simulate barbell network 

column.norm			: compute the l2 norm of a column vector 

diagSparse			: build a block diagonal matrix given block matrices
est.noise			: estimate the standard deviation of iid Gaussian noise
normalized.signal		: normalize a function to have unit standard deviation
plotCompression		: generate figure1 in the paper
one				: generate matrix of ones
row.norm			: compute the l2 norm of a row vector
sym.by.avg			: make an asymmetric matrix a symmetric one
tguhBarbell			: TGUH transformation of a noisy barbell network
uh.bu.net.inv.sm		: inverse transform with smoothing
denoise.th			: functions that denoise network signals
uh.bu.net.nonrem		: greedy version of the TGUH
uh.bu.net.nonrem.mult	: TgUH algorithm
