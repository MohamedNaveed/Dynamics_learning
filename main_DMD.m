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
    
    %control = normrnd(0,0.02);
    x(:,i+1) = model.state_propagate(i,x(:,i),control,model); % find the next state.
    
end

x_max = max(x,[],2); % max value for normalization.

%% plot the simulated data.

t_span_plot = 20;
t_step_plot = t_span_plot/model.dt + 1; 
t_idxs = (0:t_step_plot-1)*model.dt;

fig = figure;
subplot(2,1,1);
plot(t_idxs, x(1,1:t_step_plot),'LineWidth',2); 
ylabel('theta');
title('Response of the system');
subplot(2,1,2);
plot(t_idxs, x(2,1:t_step_plot),'LineWidth',2); 
ylabel('theta dot');
xlabel('time');

set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);

fig = figure;
plot(x(1,1:t_step_plot),x(2,1:t_step_plot),'LineWidth',0.5);
ylabel('theta dot');
xlabel('theta');
title('Phase portrait');


%% FFT 

fft_signal(x(:,1:t_step_plot), model, 'FFT of True data');

%% Hankel/Window DMD

window = 3; % window / time delayed samples considered for Hankel DMD.

n_samples = 60; % training samples columns of X

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
[U, S, V] = svd(X); %SVD of the data matrix 
    
 
    
thresh = 0.99999; % threshold for selecting the modes 
                   % 1 - for selecting all the modes. 
                   % 0.99999 - 99.999% energy in the modes

diag_S = diag(S);
sum_S = sum(diag_S); % sum of all the singular values. 

if thresh == 1
    A_wDMD = Xprime*pinv(X); %exact reconstruction of A (includes all the modes)
    disp('thresh = 1');
    r = length(diag_S);
else
    % finding the reduced number of modes to meet the threshold (thresh) value.
    for S_i = 1:length(diag_S)
   
        sum_S_i = sum(diag_S(1:S_i)); % partial sum of all the singular values. 
    
        if sum_S_i/sum_S >= thresh
            r = S_i;
            break;
        end      
    end
    
    A_wDMD = Xprime*V(:,1:r)*inv(S(1:r,1:r))*U(:,1:r)'; % A calculated using reduced modes.
end


fprintf('Rank of X matrix: %d \n', rank(X));
fprintf('Size of A matrix: %d \n', length(A_wDMD));
fprintf('Rank of A matrix: %d \n', rank(A_wDMD));

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

fig=figure;
subplot(2,1,1);
hold on;

plot(0:model.dt:t_span, error(1,:),'b','LineWidth',2, 'HandleVisibility','off'); 
y = ylim; % current y-axis limits
x_idx = (window)*model.dt;
plot([x_idx  x_idx],[y(1) y(2)],'k','LineWidth',2, 'DisplayName', 'Predictions start');
x_idx_train = (window + n_samples)*model.dt;
plot([x_idx_train  x_idx_train],[y(1) y(2)],'color',"#7E2F8E",'LineWidth',3, 'DisplayName', 'Training window');
ylabel('error - theta');
title('Error between wDMD and truth');
legend();

subplot(2,1,2);
hold on;
plot(0:model.dt:t_span, error(2,:),'b','LineWidth',2); 
y = ylim; % current y-axis limits
x_idx = (window)*model.dt;
plot([x_idx  x_idx],[y(1) y(2)],'k','LineWidth',2); % to show the prediction starting point
x_idx_train = (window + n_samples)*model.dt;
plot([x_idx_train  x_idx_train],[y(1) y(2)],'color',"#7E2F8E",'LineWidth',3, 'DisplayName', 'Training window');
ylabel('error - theta dot');
xlabel('time');
%saveas(fig,'error_w_20_x0_90.pdf')

fig=figure;
subplot(2,1,1);
hold on;
plot(0:model.dt:t_span, x(1,:)./x_max(1),'b','LineWidth',2); 
plot(0:model.dt:t_span, x_wDMD(1,:)./x_max(1),'--r','LineWidth',2); 
y = ylim; % current y-axis limits
x_idx = (window)*model.dt;
plot([x_idx  x_idx],[y(1) y(2)],'k','LineWidth',2);
x_idx_train = (window + n_samples)*model.dt;
plot([x_idx_train  x_idx_train],[y(1) y(2)],'color',"#7E2F8E",'LineWidth',3, 'DisplayName', 'Training window');
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
x_idx_train = (window + n_samples)*model.dt;
plot([x_idx_train  x_idx_train],[y(1) y(2)],'color',"#7E2F8E",'LineWidth',3, 'DisplayName', 'Training window');
ylabel('error - theta dot');
xlabel('time');
%saveas(fig,'wDMD_w_20_x0_90.pdf')


%% fft of error.
fft_signal(error, model, 'FFT of error');
fft_signal(x_wDMD, model, 'FFT of wDMD predictions');

 %% analyzing A

[V,D,W] = eig(A_wDMD);
diag_D = diag(D)
frequencies = logm(D(1:r,1:r))/model.dt;
diag_frequencies = diag(frequencies)./(2*pi)
plot_eigenvalues(diag_D);


