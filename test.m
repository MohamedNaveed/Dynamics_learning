%% test 
clc;clear;
N = 11;

omega_0 = 2*pi/N;

timesteps = 0:1:100;
a_vec = [1 2];
x = a_vec(1)*(exp(omega_0*1i.*timesteps) + exp(-omega_0*1i.*timesteps))...
    + a_vec(2)*(exp(2*omega_0*1i.*timesteps) + exp(-2*omega_0*1i.*timesteps));
%%
figure;
plot(timesteps,x,'LineWidth',3)

%%



