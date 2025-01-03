# Causal Discovery for Feedback Networks

Main site for the project: Sanchez-Romero, R., Ramsey, J. D., Zhang, K., Glymour, M. R., Huang, B., & Glymour, C. (2019). [Estimating feedforward and feedback effective connections from fMRI time series: Assessments of statistical methods](https://doi.org/10.1162/netn_a_00061). Network Neuroscience, 3(2), 274-306.

The paper introduced **FASK** and **Two-Step**, two new non-Gaussian causal algorithms to infer networks with feedback/cyclic interactions. The algorithms were designed and applied to fMRI data but they can be used for other non-Gaussian, linear datasets. See paper for details.

## Non-Gaussian Algorithms

Fast Adjacency Skewness (**FASK**): available as part of the Tetrad project, https://github.com/cmu-phil/tetrad

**Two-Step**: available as a standalone Matlab repository, https://github.com/cabal-cmu/Two-Step 

## Datasets

Simulated and empirical data used in the paper. 

Download here https://bit.ly/datafeedbacks

The data folder contains various sub-folders:

1. Simple Networks: Synthetic data for 9 cyclic networks (Figure 1 in the text) and their amplifying and control variants, resulting in 18 simulations. For each simulation: 60 individual datasets and 60 datasets of ten concatenated invididuals. 

2. Macaque Networks: Synthetic data for 4 cyclic networks based on Macaque structural connections from Markov et al.,(2013). For each simulation: 60 individual datasets and 60 datasets of ten concatenated invididuals.

3. Measurement Noise Data: Simulated data with and without measurement noise for the analysis of Section 4.3.1., and Figure 4. 

4. Temporal Undersampling Data: Simulated data sampled at different resolutions for the analysis of Section 4.3.2., and Figure 5.

5. Resting-State fMRI Data from the Medial Temporal Lobe (MTL): Empirical resting-state data from Shah et al.,(2017). We include here the data as we processed it for our analysis as described in the text. We include 23 individual datasets and 23 datasets of ten concatenated individuals. The original data is available at https://github.com/shahpreya/MTLnet.

6. Rhyming Task fMRI Data: Empirical task data preprocessed from Ramsey et al.,(2010). We include 9 individual datasets and 1 dataset of nine concatenated individuals. The original raw NIfTI files can be downloaded at https://openfmri.org/dataset/ds000003/

7. Parameter Tuning Data for Simple Networks: Synthetic data for Networks 5 and 7 for parameter tuning of the algorithms as described in Supplementary Material B.
