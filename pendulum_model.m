function [model] = pendulum_model()

model.m = 0.5; % mass
model.dt = 0.1; % discretization time-step
model.nx = 2; % state dimension
model.nu = 1; % control dimension
model.g = 9.81; % gravity constant
model.L = model.g/(4*(pi^2)); % length. chosen specifically for the time period to be 1 for small angles.  
model.state_propagate = @pendulum_nl_state_prop; % state propagation function/dynamics.
model.c = 0; %damping coefficient. 
end