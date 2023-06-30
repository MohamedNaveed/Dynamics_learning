% main file 
clc;clear;

%% define pendulum
model = pendulum_model();
x0 = [deg2rad(90),0]; %initial state. 0 - pendulum hanging downwards. 

%% simulate the pendulum motion

t_span = 100; % (s)
t_steps = t_span/model.dt;
x = zeros(model.nx, t_steps+1);
x(:,1) = x0;
control = 0;

for i = 1:t_steps
    
    x(:,i+1) = pendulum_nl_state_prop(i,x(:,i),control,model);
    
end

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

A = Xprime*pinv(X);


%% prediction using A

x_DMD = zeros(model.nx,t_steps+1);
x_DMD(:,1) = x0;

for i=1:t_steps

    x_DMD(:,i+1) = A*x_DMD(:,i); 

end

%% error 

error = x - x_DMD;
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
plot(0:model.dt:t_span, x(1,:),'b','LineWidth',2); 
plot(0:model.dt:t_span, x_DMD(1,:),'--r','LineWidth',2); 
ylabel('error - theta');
legend('True','DMD');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, x(2,:), 'b','LineWidth',2); 
plot(0:model.dt:t_span, x_DMD(2,:),'--r','LineWidth',2);
ylabel('error - theta dot');
xlabel('time');

%% window DMD

window = 6;
n_samples = 500;

X = zeros(model.nx*window,n_samples);

for w = 1:window
    
    X(model.nx*(w-1) +1:model.nx*w,:) = x(:,window - (w - 1) : window - (w - 1) + n_samples -1);
    
end

Xprime = zeros(model.nx*window,n_samples);

for w = 1:window
    
    Xprime(model.nx*(w-1) +1:model.nx*w,:) = x(:,window - (w - 1) + 1 : window - (w - 1) + n_samples );
    
end

% solve for A 

A = Xprime*pinv(X);
disp('rank of X')
rank(X)
disp('length of A')
length(A)
%% prediction using A

x_wDMD = zeros(model.nx,t_steps+1);

y_wDMD = zeros(model.nx*window,1); %observable

%initial condition
for w = 1:window
    y_wDMD(model.nx*(window - w) + 1:model.nx*(window - w + 1)) = x(:,w);
end

x_wDMD(:,1:window) = x(:,1:window); %assume same for the first few steps. 

for i = window:t_steps

    y_wDMD = A*y_wDMD; 
    x_wDMD(:,i+1) = y_wDMD(1:model.nx);
end

%% error 

error = x - x_wDMD;
%% plot the data.

figure(4);
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, error(1,:),'b','LineWidth',2); 
ylabel('error - theta');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, error(2,:),'b','LineWidth',2); 
ylabel('error - theta dot');
xlabel('time');

figure(5);
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, x(1,:),'b','LineWidth',2); 
plot(0:model.dt:t_span, x_wDMD(1,:),'--r','LineWidth',2); 
ylabel('error - theta');
legend('True','wDMD');
subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, x(2,:), 'b','LineWidth',2); 
plot(0:model.dt:t_span, x_wDMD(2,:),'--r','LineWidth',2);
ylabel('error - theta dot');
xlabel('time');

