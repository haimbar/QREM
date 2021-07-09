# QREM - Quantile Regression EM algorithm

The estimating equations for quantile regression (QR) can be solved using an EM algorithm in which the M-step is computed via weighted least squares, with weights computed at the E-step as the expectation of independent generalized inverse-Gaussian variables.
The package extends QR to allow for random effects in the linear predictor. The generalized alternating minimization (GAM) framework also allows for variable selection in the QR, large P setting. 

For more details on the theory, see the paper on [ArXiV](https://arxiv.org/abs/2104.08595). For more usage documentation, see the [QREM.pdf](doc/QREM.pdf) file in the doc/ subfolder in this repository .
