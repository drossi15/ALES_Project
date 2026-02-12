# ALES_Project
 Unibg Adaptive Learning, Estimation and Supervision of Dynamical Systems course project

This repository presents three different strategies for Recursive Principal Component Analysis (RPCA), based on the paper:

Li et al., Recursive PCA for Adaptive Process Monitoring
https://www.sciencedirect.com/science/article/pii/S0959152400000226

Implemented Approaches

The project includes the following RPCA update strategies:

1️⃣ Full Update (Sample-by-Sample)

Updates the full correlation structure at each new sample

Uses the VRE (Variance of Reconstruction Error) criterion to select the number of principal components

2️⃣ Rank-One Update (Sample-by-Sample)

Performs a rank-one update for each incoming sample

Uses the CPV (Cumulative Percentage of Variance) criterion to select the principal components

3️⃣ Block Update + Lanczos Method

Processes data in blocks instead of single samples

Uses the Lanczos algorithm for efficient eigendecomposition

Applies the CPV criterion for principal component selection
