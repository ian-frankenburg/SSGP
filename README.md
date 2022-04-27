# spSSGP Example
To run an example analysis using the nonlinear partial differential equation provided in the manuscript, first compile the ```spatialCalibrator.cpp``` function housed within the ```cpp models``` folder. Once this is compiled, the environment will have access to the necessary function call to fit the model.

The ```helper.R``` file will generate the synthetic PDE data used for this demo, but the only code to be run is the ```demo.R``` file. Once the model is fit, prior-to-posterior learning plots are generated within the same file.

This example illustrates nonlinear PDE inference on a 10 x 10 grid for simplicity. We make use of 25 training runs where the numerical solution to the PDEs are computed across a Latin hypercube design space. The MCMC chain is run for 1000 iterations by setting ```niter=1000``` with the first 500 discarded for warm-up, illustrating the fast convergence of the FFBS algorithm.
