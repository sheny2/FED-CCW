# G-computation oracle formulation

Let site membership be \(K\in\{1,\ldots,5\}\), sampled with

\[
P(K=k)=\frac{n_k}{\sum_j n_j}.
\]

Conditional on site, baseline covariates are generated from the configured
site-specific distributions. At interval \(m\), the longitudinal covariates
follow

\[
L_{jm}=\eta_{0j}+\eta_{1j}L_{j,m-1}+\eta_{2j}X_1+
\eta_{3j}X_2+\epsilon_{jm},\qquad j\in\{1,2\}.
\]

For the initiate-within-grace-period strategy, treatment follows the true
conditional initiation hazard before the last grace-period interval:

\[
P(A_m=1\mid A_{m-1}=0,X,L_m,K=k)
=\operatorname{expit}\{\alpha_k+\gamma^T(X,L_m)\},\qquad m<\tau.
\]

If treatment has not started previously, \(A_\tau=1\) is imposed. Treatment
then remains on. Under the no-initiation strategy, \(A_m=0\) throughout.

For strategy \(g\), the true conditional event hazard is

\[
p_m^g=\operatorname{expit}\{\beta_0+\beta_X^TX+
\beta_L^TL_m+\beta_{L^-}^TL_{m-1}+\beta_{\mathrm{trt}}A_m^g\}.
\]

For each Monte Carlo trajectory, conditional survival is evaluated as

\[
S_i^g(m)=\prod_{j=1}^m(1-p_{ij}^g),
\]

and the standardized survival curve is

\[
S^g(m)=E\{S_i^g(m)\},
\]

where the expectation is over the scenario-specific site mixture, baseline
covariates, longitudinal covariates, and stochastic pre-tau initiation timing.
The oracle reports risks, RD, RR, OR, RMST in each arm, and the RMST
difference from these standardized survival curves.

