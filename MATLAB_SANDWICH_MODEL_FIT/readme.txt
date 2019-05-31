fitme is the script that starts the whole process of fitting (although it assumed that the real input-output data + confidence intervals where available in the workspace)
It calls iteratively the function 'anneal()', which is my implementation of a montecarlo optimization (i.e. simulated annealing - see Kirkpatrick et al. seminal work; or refer to Numerical Recipes).

Anneal() calls cost() to determine whether the actual change of the (free) model parameters to fit lead to an improvement in the fitness index..

cost() is the actual fitness index that is calculated at any step of the optimization
It calls 'Model()'..

Model() is the actual routine that performs the model 'simulation' (which is actually frequency-domain filtering)..

plot_model   is a script to simply provide a visual indication of the fit performance (or to plot the model output)
