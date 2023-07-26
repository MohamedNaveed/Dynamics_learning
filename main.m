% main file 
clc;clear;

%% define pendulum      
model = pendulum_model(); 
x0 = [deg2rad(90),0]; %initial state (deg) . 0 - pendulum hanging downwards. 

%% simulate the pendulum motion (ground truth data)

t_span = 20; % (seconds) 
t_steps = t_span/model.dt; % number of time steps
x = zeros(model.nx, t_steps+1); % state

x(:,1) = x0;
control = 0; % No control actions. 

for i = 1:t_steps
    
    x(:,i+1) = model.state_propagate(i,x(:,i),control,model); % find the next state.
    
end


x_max = max(x,[],2); % max value for normalization.


%% plot the simulated data.

fig = figure(1);
subplot(2,1,1);
plot(0:model.dt:t_span, x(1,:),'LineWidth',2); 
ylabel('theta');
title('Response of the system');
subplot(2,1,2);
plot(0:model.dt:t_span, x(2,:),'LineWidth',2); 
ylabel('theta dot');
xlabel('time');

set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);

%% exact DMD

X = x(:,1:end-1);  % data matrix

Xprime = x(:,2:end);
 
% solve for A 

A_DMD = Xprime*pinv(X); % finding the A matrix that fits the data. (keeping all the modes)


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
title('Error between DMD predictions and true data');

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
title('Predicted response using DMD');

subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, x(2,:)./x_max(2), 'b','LineWidth',2); 
plot(0:model.dt:t_span, x_DMD(2,:)./x_max(2),'--r','LineWidth',2);
ylabel('error - theta dot');
xlabel('time');


%% Hankel/Window DMD

window = 20; % window / time delayed samples considered for Hankel DMD.
n_samples = 81; % training samples columns of X

X = zeros(model.nx*window,n_samples);

% creating the data matrix.
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
fprintf('Rank of X matrix: %d \n', rank(X));
fprintf('Size of A matrix: %d \n', length(A_wDMD));

%% calculate error in the training data.

error_fit = Xprime - A_wDMD*X;

fprintf('Max training error using A = %d \n', max(max(error_fit)));

%% prediction using A

x_wDMD = zeros(model.nx,t_steps+1);

y_wDMD = zeros(model.nx*window,1); % observable. stacked observable for Hankel DMD.

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

fig=figure(4);
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, error(1,:),'b','LineWidth',2); 
y = ylim; % current y-axis limits
x_idx = (window)*model.dt;
plot([x_idx  x_idx],[y(1) y(2)],'k','LineWidth',2);
ylabel('error - theta');
title('Error between wDMD predictions and true data');

subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, error(2,:),'b','LineWidth',2); 
y = ylim; % current y-axis limits
x_idx = (window)*model.dt;
plot([x_idx  x_idx],[y(1) y(2)],'k','LineWidth',2); % to show the prediction starting point
ylabel('error - theta dot');
xlabel('time');
%saveas(fig,'error_w_20_x0_90.pdf')

fig=figure(5);
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, x(1,:)./x_max(1),'b','LineWidth',2); 
plot(0:model.dt:t_span, x_wDMD(1,:)./x_max(1),'--r','LineWidth',2); 
y = ylim; % current y-axis limits
x_idx = (window)*model.dt;
plot([x_idx  x_idx],[y(1) y(2)],'k','LineWidth',2);
ylabel('error - theta');
legend('True','wDMD');
title('Predicted response using wDMD');

subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, x(2,:)./x_max(2), 'b','LineWidth',2); 
plot(0:model.dt:t_span, x_wDMD(2,:)./x_max(2),'--r','LineWidth',2);
ylim([-1,1])
y = ylim; % current y-axis limits
x_idx = (window)*model.dt;
plot([x_idx  x_idx],[y(1) y(2)],'k','LineWidth',2);
ylabel('error - theta dot');
xlabel('time');
%saveas(fig,'wDMD_w_20_x0_90.pdf')

%% Comparison with Autoregressive model.
%{ 
 

%% Autoregressive AR model. 

X = zeros(model.nx*window,n_samples);

% creating data matrix
for w = 1:window
    
    X(model.nx*(w-1) +1:model.nx*w,1:n_samples) = x(:,window - (w - 1) : ...
                                           window - (w - 1) + n_samples -1);
end

Xprime_arma = x(:,window + 1 : window + n_samples );
    
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
title('Error between AR predictions and true data');

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
%}