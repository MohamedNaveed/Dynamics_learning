function [model] = pendulum_model()

model.m = 0.5;
model.u_max = 1;
model.dt = 0.1;
model.nx = 2;
model.nu = 1;
model.g = 9.81;
model.L = model.g/(4*(pi^2)); % for the time period to be 1 for small angles.  
end