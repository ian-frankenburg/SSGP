# spSSGP Example
To run an example analysis using the nonlinear partial differential equation provided in the manuscript, first compile the ```spatialCalibrator.cpp``` function housed within the ```cpp models``` folder. Once this is compiled, the environment will have access to the necessary functions.

The ```helper.R``` file will generate the synthetic PDE data used for this demo, but the only code to be run is the ```demo.R``` file. Once the model is fit, a prior-to-posterior learning plots are generated.

This example illustrates nonlinear PDE inference on a 10 x 10 grid for simplicity. We make use of 25 training runs where the numerical solution to the PDEs are computed across a Latin hypercube design space. The MCMC chain is run for 1000 iterations with the first 500 discarded for warm-up, illustrating the fast convergence of the FFBS algorithm.
