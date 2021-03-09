1. "sim_main" is the main file used to run the demo code. It consists of five main steps:
   (1) Set simulation parameters;
   (2) Generate a dataset;
   (3) Fit the model;
   (4) Policy search;
   (5) Store the output.

2. "function" defines functions that are needed in above procedures.
   (1)  "corfx": computes the values of the Matern correlation function
	Arguments: 
	  d: distance matrix;
	  theta: parameters for Matern correlation.
	Value produced:
	  Matern correlation matrix.

   (2)  "ST_Krige": the function of fitting the Spatiotemporal model.
	Arguments:
	  Y: disease rate (logit transform);
	  X: environmental covariates;
	  trt: resource allocation;
	  A: adjancy matirx;
	  a, b, sd_beta: hyperparameters;
	  iters: number of iterations in the MCMC sampling;
	  burn: number of burn-in samples.
	Value produced:
	  the posterior samples of the parameters.

   (3)  "policy_highY": policy that allocate all resources to healthzones with highest prevalence;
	"policy_even": policy that allcate all resources evenly to all health zones;
	"policy_no: policy that allocate no resources;
	"policy_all": policy that allocate resources to everyone;
	"policy_linear":policy that we proposed with linear utility function;
	"policy_quad": policy that we proposed with quadratic utility function.	
	Arguments: 
	  X: environmental covariates;
	  Y: disease rate;
	  nb_Y: average disease rate of neighbor zones;
 	  alpha: risk factor weights.
	Value produced: recommended resource allocation;

   (4)  "get.value.based.pos": estimates the loss value based on fitted model.
	Arguments:
	  alpha: risk factor weights;
	  params: parameters of the model;
	  X_env: environmental covariates;
	  Y_base: baseline disease rate (logit transform);
	  nyear: number of years to simulate forward to calculte the loss value;
	  nsamps: number of samples to simulate;
	  policy: the policy used to allocate resources, it can be "policy_linear", "policy_quad", etc.
	Value produced: The estimated loss value associated with the policy.

   (5)  "get.value.true": estimates the loss value based on the true model.
	Arguments:
	  alpha: risk factor weights;
	  params: parameters of the model;
	  X_env: environmental covariates;
	  Y_base: baseline disease rate (logit transform);
	  latent_base: the spatiotemporal process baseline value;
	  nyear: number of years to simulate forward to calculte the loss value;
	  nsamps: number of samples to simulate;
	  policy: the policy used to allocate resources, it can be "policy_linear", "policy_quad", etc.
	Value produced: The true loss value associated with the policy.

 3. "DRC_malaria_data_analysis" contains the code for DRC malaria data analysis. It includes reading related data,
    preprocessing the data, fitting the spatiotemporal model and estimating one optimal policy averaging over 
    the uncertainty of parameters.

 4. "DRC_malaria_data_analysis_posterior" contains the code for estimating posterior distribution of optimal policy parameters
    for DRC malaria data analysis. It includes reading related data, preprocessing the data, fitting the spatiotemporal model,
    and estimating the posterior distribution of the optimal policy parameters (alpha_opt).