# Data Driven modeling using Dynamic Mode Decomposition (DMD).

The repository has implementations of the basic/exact DMD and Hankel DMD for the nonlinear pendulum model. 

main.m implements DMD for modeling the trajectory originating from a single initial condition.

main_ensemble_training.m implements DMD for modeling a set (ensemble) of trajectories originating from initial conditions from a domain.

A different system can be used by defining the model as done for the pendulum system and calling that function. 
