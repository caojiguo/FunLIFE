# The R codes for simulation studies and real data analysis in the manuscript entitled "Multi-Resolutional Panel Data Models with Interactive Fixed Effects". 

###########################################################

1. Funlife_est.R: the R code used to calculate the estimation for the parameters of the functional linear regression with iterative fixed effects(FunLife). The initial output is used as the starting point of the iteration algorithm for computing the estimations of the proposed method.

2. Funlife_fact.R: the R code used to identify the optimal factor for the proposed model when the true factor number is not given. 

3. Funlife_para.R: the R code used to choose the smoothing parameter for the penalty for functional coefficients in the proposed model.

4. FLM_est.R: the R code used to calculate the conventional functional linear model without considering the iterative fixed effects in the proposed model.

5. FLM_eval.R: the R code used to evaluate the the performance of estimates in both proposed Funlife model and conventional functional linear mode without considering the iterative fixed effects.

6. Pdm_est.R: the R code used to calculate the estimates of conventional panel data model using the same setting in the proposed model.

7. Sam_data.R: the R code used to generate data for the simulation studies presented in the manuscript.

###########################################################
