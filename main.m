% main file 
clc;clear;

%% define pendulum  
model = pendulum_model(); 
x0 = [deg2rad(90),0]; %initial state. 0 - pendulum hanging downwards. 15

%% simulate the pendulum motion

t_span = 10; % (s)
t_steps = t_span/model.dt;
x = zeros(model.nx, t_steps+1);

%sample path 1
x(:,1) = x0;
control = 0;

for i = 1:t_steps
    
    x(:,i+1) = pendulum_nl_state_prop(i,x(:,i),control,model);
    
end

%sample path 2
x0 = [deg2rad(60),0];
x_2 = zeros(model.nx, t_steps+1);
x_2(:,1) = x0;
control = 0;

for i = 1:t_steps
    
    x_2(:,i+1) = pendulum_nl_state_prop(i,x_2(:,i),control,model);
    
end



x_max = max(x,[],2); % max value for normalization.

%{
%% plot the data.

figure(1);
subplot(2,1,1);
plot(0:model.dt:t_span, x(1,:),'LineWidth',2); 
ylabel('theta');
subplot(2,1,2);
plot(0:model.dt:t_span, x(2,:),'LineWidth',2); 
ylabel('theta dot');
xlabel('time');

%% simple DMD

%data matrix


X = x(:,1:end-1); 
Xprime = x(:,2:end);
 
% solve for A 

A_DMD = Xprime*pinv(X);


%% prediction using A

x_DMD = zeros(model.nx,t_steps+1);
x_DMD(:,1) = x0;

for i=1:t_steps

    x_DMD(:,i+1) = A_DMD*x_DMD(:,i); 

end

%% error 

error = (x - x_DMD)./x_max; %normalized

%% plot the data.

figure(2);
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, error(1,:),'b','LineWidth',2); 
ylabel('error - theta');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, error(2,:),'b','LineWidth',2); 
ylabel('error - theta dot');
xlabel('time');

figure(3);
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, x(1,:)./x_max(1),'b','LineWidth',2); 
plot(0:model.dt:t_span, x_DMD(1,:)./x_max(1),'--r','LineWidth',2); 
ylabel('error - theta');
legend('True','DMD');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, x(2,:)./x_max(2), 'b','LineWidth',2); 
plot(0:model.dt:t_span, x_DMD(2,:)./x_max(2),'--r','LineWidth',2);
ylabel('error - theta dot');
xlabel('time');

%% monte_carlo test

[error_mean_DMD, error_std_DMD] = monte_carlo_test(A_DMD, 'DMD', model);

figure;
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, error_mean_DMD(1,:),'b','LineWidth',2);
plot(0:model.dt:t_span, error_mean_DMD(1,:) + error_std_DMD(1,:),'--r','LineWidth',2);
plot(0:model.dt:t_span, error_mean_DMD(1,:) - error_std_DMD(1,:),'--r','LineWidth',2);
ylabel('error - theta');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, error_mean_DMD(2,:),'b','LineWidth',2);
plot(0:model.dt:t_span, error_mean_DMD(2,:) + error_std_DMD(1,:),'--r','LineWidth',2);
plot(0:model.dt:t_span, error_mean_DMD(2,:) - error_std_DMD(1,:),'--r','LineWidth',2);
 
ylabel('error - theta dot');
xlabel('time');
%}

%% window DMD

window = 10;
n_samples = 50;

X = zeros(model.nx*window,n_samples);

for w = 1:window
    
    X(model.nx*(w-1) +1:model.nx*w,:) = x(:,window - (w - 1) : ...
                                           window - (w - 1) + n_samples -1);
    
end

Xprime = zeros(model.nx*window,n_samples);

for w = 1:window
    
    Xprime(model.nx*(w-1) +1:model.nx*w,:) = x(:,window - (w - 1) + 1 : ...
                                            window - (w - 1) + n_samples );
    
end

% solve for A 

A_wDMD = Xprime*pinv(X);
disp('rank of X')
rank(X)
disp('length of A')
length(A_wDMD)

%% prediction using A

x_wDMD = zeros(model.nx,t_steps+1);

y_wDMD = zeros(model.nx*window,1); %observable

%initial condition
for w = 1:window
    y_wDMD(model.nx*(window - w) + 1:model.nx*(window - w + 1)) = x(:,w);
end

x_wDMD(:,1:window) = x(:,1:window); %assume same for the first few steps. 

for i = window:t_steps

    y_wDMD = A_wDMD*y_wDMD; 
    x_wDMD(:,i+1) = y_wDMD(1:model.nx);
end

%% error 

error = (x - x_wDMD)./x_max;
%% plot the data.

figure(4);
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, error(1,:)./x_max(1),'b','LineWidth',2); 
ylabel('error - theta');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, error(2,:)./x_max(2),'b','LineWidth',2); 
ylabel('error - theta dot');
xlabel('time');

figure(5);
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, x(1,:)./x_max(1),'b','LineWidth',2); 
plot(0:model.dt:t_span, x_wDMD(1,:)./x_max(1),'--r','LineWidth',2); 
ylabel('error - theta');
legend('True','wDMD');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, x(2,:)./x_max(2), 'b','LineWidth',2); 
plot(0:model.dt:t_span, x_wDMD(2,:)./x_max(2),'--r','LineWidth',2);
ylabel('error - theta dot');
xlabel('time');

%% ARMA model. 1 sample

X = zeros(model.nx*window,2*n_samples);

for w = 1:window
    
    X(model.nx*(w-1) +1:model.nx*w,1:n_samples) = x(:,window - (w - 1) : ...
                                           window - (w - 1) + n_samples -1);
    X(model.nx*(w-1) +1:model.nx*w,n_samples+ 1: 2*n_samples) = x_2(:,window - (w - 1) : ...
                                           window - (w - 1) + n_samples -1);

end

Xprime_arma = [x(:,window + 1 : window + n_samples ), ...
                x_2(:,window + 1 : window + n_samples )];
    
A_arma = Xprime_arma*pinv(X);

%% prediction using A

x_arma = zeros(model.nx,t_steps+1);

y_arma = zeros(model.nx*window,1); %observable

%initial condition
for w = 1:window
    y_arma(model.nx*(window - w) + 1:model.nx*(window - w + 1)) = x(:,w);
end

x_arma(:,1:window) = x(:,1:window); %assume same for the first few steps. 

for i = window:t_steps

    x_arma(:,i+1) = A_arma*y_arma; 
    y_arma = [x_arma(:,i+1);y_arma(1:model.nx*(window-1))];
end

%% error 

error = (x - x_arma)./x_max;
%% plot the data.

figure(6);

subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, error(1,:)./x_max(1),'b','LineWidth',2); 
ylabel('error - theta');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, error(2,:)./x_max(2),'b','LineWidth',2); 
ylabel('error - theta dot');
xlabel('time');

%% error 

error = (x_wDMD - x_arma)./x_max;
%% plot the data.

figure(7);

subplot(2,1,1);
title('Error between ARMA and wDMD');
hold on;
plot(0:model.dt:t_span, error(1,:)./x_max(1),'b','LineWidth',2); 
ylabel('error - theta');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, error(2,:)./x_max(2),'b','LineWidth',2); 
ylabel('error - theta dot');
xlabel('time');

%% ensemble data fit - training
n_samples = 100;
window = 10;
[A_arma, error_fit_arma] = model_fit('ARMA', model, window, n_samples);

%% monte_carlo test - testing
n_mc_runs = 100;
[error_mean_arma, error_std_arma] = monte_carlo_test(A_arma, 'ARMA', model,...
                                    window, n_mc_runs);

%%
figure;
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, error_mean_arma(1,:),'b','LineWidth',2);
plot(0:model.dt:t_span, error_mean_arma(1,:) + error_std_arma(1,:),'--r','LineWidth',2);
plot(0:model.dt:t_span, error_mean_arma(1,:) - error_std_arma(1,:),'--r','LineWidth',2);
ylabel('error - theta');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, error_mean_arma(2,:),'b','LineWidth',2);
plot(0:model.dt:t_span, error_mean_arma(2,:) + error_std_arma(1,:),'--r','LineWidth',2);
plot(0:model.dt:t_span, error_mean_arma(2,:) - error_std_arma(1,:),'--r','LineWidth',2);
 
ylabel('error - theta dot');
xlabel('time');

