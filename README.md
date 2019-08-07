
# Modeling Syphilis in Louisiana and Massachusetts 

![modeled natural history](inst/figures/background_figures/modeled natural history.png?raw=true)

This is a deterministic compartmental model of syphilis implemented in R using Rcpp. 

The project workflow is as follows: 

### Build the model 

This step has already been done for a user of the package. 
If you need to make changes to the model, this is where you would start. 

To build the model, one must setup structural equations, prior distributions, data
targets, and likelihood distributions. 

These are defined (mostly) in: 

  - `src/syph_trans_model.cpp` (differential equations), 

  - `inst/extdata/priors.txt` (prior distribution specifications),

  - `R/model.priors.fun.R` (prior distribution density function), 
    data target files in `inst/extdata/` (read in with `R/initial.files.R`), 

  - `R/model.likelihood.fun.R` (the likelihood distribution density constructed from model simulation 
    which measures the likelihood of data targets given the hypothesized model). 

To load the model into R, run the following: 

    library(syphilis) # after using `devtools::install()` to install the package
    state <- 'LA' # or 'MA'
    load.start() # load the model parametrization

### Calibrate the model

We use a Bayesian approach to estimating parametric
uncertainty associated with fitting our model to the data targets. We first
optimize the model to get best model fit, and then start Monte Carlo Markov
Chains from parameter vectors which are the outcomes of several model
optimizations. 

We use an adaptive MCMC algorithm to better approximate the posterior
distribution given that we expect many of the model parameters to be
correlated (or inversely correlated) with each other in the posterior
distribution.

The code to optimize the model is in `analysis/optimize/`. We chose to use a
simultaneous calibration approach where certain parameters (those controlling
the parameterization of the natural history of syphilis) were shared for the
Louisiana and Massachusetts model parametrizations. This approach is 
implemented in
`analysis/optimize/construct_simultaneous_optimization_function.R`.

The optimization code is designed to be run on a cluster using the
[SLURM](https://slurm.schedmd.com/) workload scheduler. See the sbatch 
scripts here: `analysis/optimize/run_script.sh`. Notice that the optimization 
is run in an array of 25 to yield 25 model optimizations. 

Part of the model building process involves confirming (qualitatively) that the 
optimizations land in sufficiently similar places both with respect to model 
parameters and the log-posterior values (computed by the `dLogPosterior`) 
function. This is a necessary step to establish that the posterior distribution 
is shaped such that the global optimum (best model fit) can be reliably found 
and indicates that the MCMC process for sampling parametric uncertainty may 
be able to perform well (though this is not necessarily guaranteed). 

The code to run the MCMC sampling algorithm is here:
`analysis/mcmc/run_amcmc_from_optims.R`.  This code is also designed to be
scheduled using slurm. See its sbatch scripts in the same folder.

Notice that the optimizations and MCMC sampling code store their outputs 
in the `inst/optims/` and `inst/mcmc/` folders. This allows code inside the
syphilis project may access those files by using the `system.file(...,
package='syphilis')` syntax.

### Simulate Intervention Scenarios
 
In order to simulate hypothetical future interventions, we change the 
future values of the model estimated screening rates from the assumption 
of a continuation of screening at the 2016 level. 

This is implemented in `R/screening_interventions.R`.



