# CACE
The method is to estimate the mean and site-specific CACE from a multisite randomized trial with binary outcomes and one-sided noncompliance (Little and Yau 1998) by estimating a bivariate random coefficient model via maximum likelihood (ML).

There are four demonstration scripts (main scripts) uploaded to the repository. They are r1_q4_n40_x1x2_potential.R, r1_q4_n80_J65_x1x2_my05.R, r1_q4_eassist_x1x2_temp.R and r3_q4_n40_J170_x1x2_my0.R. The rest of the .R files in the repository are R function scripts. Users need to download all the R function scripts to run any of the main scripts. The main script contains what R functions required to be imported into the R environment. Script r1_q4_n40_x1x2_potential.R simulates a dataset from a multisite study with cluster size of 40 nested within 170 sites, with potential outcomes and observed compliances, and estimates the average treatment effect by lme4 package. This is to show the gold standard estimation of CACE with potential outcomes. The simulation function in the script generates a dataset with potential or observed outcomes as explained by the comments in the R script. The datasets of all potential outcomes and observed outcomes are named as "L2" and "L1o" respectively. Script r1_q4_n80_J65_x1x2_my05.R simulates a dataset from a multisite study with cluster size of 80 nested within 65 sites and 5% missing outcomes, and estimates the CACE by the method with never taker's random effect shared by others. Script r1_q4_eassist_x1x2_temp.R simulates a dataset from the e-assist study with an average cluster sample size of 11 nested within 170 physicians. The treatment effect is estimated with never taker's random effect shared with others and 4 abscissas for AGHQ. Script r3_q4_n40_J170_x1x2_my0.R simulates a dataset from a multisite study with a fixed cluster sample size of 40 nested within 170 sites and completely observed outcomes, and estimates the CACE by the method with no shared random effect. Users can modify the main script according to their datasets and needs, such as number of abscissas for AGHQ (Q), reduced dimension (r_prime), total number of iterations (niter) etc. 

The function in script cace_1side_parallel_lambda.R is the major function for CACE estimation with shared random effects (r=1 or 2). The function in script cace_1side_parallel_lambda.R is the major function for CACE estimation with shared random effects (r=1 or 2).The function in script cace_1side_parallel.R is the major function for CACE estimation without shared random effects (r=3). When use the function, the user needs to define J (number of sites), but doesn’t need to define n (#subjects nested in a site) because the function calculates n for each site. For one-sided noncompliance, r should be 3 and k be 1, and users just set r=3 and k=1 for all one-sided noncompliance cases. r_prime is the number of random effect remained after random effect sharing. I usually set r_prime=1 to save time. In most cases, r_prime=1 should work and give results. L1o in the main script is the dataset name. x0 and x1 are the vectors of covariates controlled in Y and C models, and they are set as characters. init is where to put initial values. col.clinic, col,trt, col.D and col.Y are the column numbers of variables cluster/site, treatment assignment, treatment receipt and outcome, and they are numeric.

The function in script set_init.R generates the initial values of all parameters when there is one or two shared random effect(s) (r=1 or 2), while the function in script set_init_noshare.R generates the initial values of all parameters when there is no random effect sharing (r=3). The functions in the rest of the R scripts in function folder can be used with or without shared random effect (r=1, 2 or 3). For set_init(), when it is one-sided noncompliance, set r=3 and side=1. In the main script, C ~ x1+x2+(1|clinic), Y ~ x1+x2+(1|clinic) are the formulas for the logsitic model fitting. “clinic” is the name of the cluster level variable. “C” and “Y” are the variable names of compliance and outcome. "ID" is the subject level variable name and "pmm" means predictive mean matching. r_prime is the number of random effects after random effect sharing, and it is set to one mostly. col.clinic, col,trt, col.D and col.Y are the column numbers of variables cluster/site, treatment assignment, treatment receipt and outcome, and they are numeric. Please keep initial values in this order: alpha, gamma, lambda, tau, delta. alpha is the fixed effect, lambda is the coefficient for random effect sharing, tau is covariance matrix elements in Y model. gamma is the fixed effect and delta is variance in C model.

For successful running of functions cace_1side_parallel() or cace_1side_parallel_lambda(), the datasets must contain cluster ID, treatment assignment (numeric, 0=assigned to control or 1=assigned to treatment), treatment receipt (numeric, 0=receiving control or 1=receiving treatment) and outcome variables (numeric, 0 or 1). The cluster ID, treatment assignment and receipt variables must be completely observed. The outcome variable may be completely observed or missing at random. The function can handle the missing outcomes MAR.
