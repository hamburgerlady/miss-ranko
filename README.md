# miss-ranko
Robust low rank matrix factorization with missing data.
Matlab implementation of RANSAC framework based on minimal solvers 
for low rank factorization of matrices with missing data.

Please cite 
M. Oskarsson, K. Batstone and K. Åström "Trust No One: Low Rank Matrix Factorization Using Hierarchical RANSAC"
in Proc. CVPR 2016, Las Vegas, USA
F. Jiang, M. Oskarsson and K. Åström "On the minimal problems of low-rank matrix factorization" 
in Proc. CVPR 2015, Boston, USA

Contains two versions
* genrank-solver - Solver for a general fixed rank
* sfmrank4-solver - Solver for affine structure from motion, using rank 4 factorization