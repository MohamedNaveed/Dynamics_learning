% modeling periodic signals using AR model.
clc;clear;

%% define signal
fs = 100; % sampling freq.
dt = 1/fs; % sampling interval
model.dt = dt;
model.nx = 1;
t_span = 10;
t = 0:dt:t_span;
t_steps = t_span/dt;
f1 = 1;
f2 = 10;

x = sin(2*pi*f1*t) + 0.5*sin(2*pi*f2*t);
%x = sin(2*pi*f1*t);
x_max = max(x);
figure;
plot(t, x, 'LineWidth',2);
ylabel('signal');
xlabel('time');
title('true signal');

%% plot fft.

fft_signal(x,model,'FFT of True signal');

%% fit AR model.

window = 4; % window / time delayed samples considered.
n_samples = 81; % training samples columns of X

X = zeros(model.nx*window,n_samples);

% creating the data matrix.
for w = 1:window
    
    X(model.nx*(w-1) +1:model.nx*w,:) = x(:,window - (w - 1) : ...
                                           window - (w - 1) + n_samples -1);
    
end

[U, S, V] = svd(X); %SVD of the data matrix 

Xprime_arma = x(:,window + 1 : window + n_samples );

thresh = 0.99999; % threshold for selecting the modes 
                   % 1 - for selecting all the modes. 
                   % 0.99999 - 99.999% energy in the modes
if size(S,1) == 1
    diag_S = S(1);
else
    diag_S = diag(S);
end
sum_S = sum(diag_S); % sum of all the singular values. 

if thresh == 1
    A_arma = Xprime_arma*pinv(X); %exact reconstruction of A (includes all the modes)
else
    % finding the reduced number of modes to meet the threshold (thresh) value.
    for S_i = 1:length(diag_S)
   
        sum_S_i = sum(diag_S(1:S_i)); % partial sum of all the singular values. 
    
        if sum_S_i/sum_S >= thresh
            r = S_i;
            break;
        end      
    end
    
    A_arma = Xprime_arma*V(:,1:r)*inv(S(1:r,1:r))*U(:,1:r)'; % A calculated using reduced modes.
end

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

figure;
hold on;
plot(0:model.dt:t_span, error(1,:),'b','LineWidth',2, 'HandleVisibility','off'); 
y = ylim; % current y-axis limits
x_idx = (window)*model.dt;
plot([x_idx  x_idx],[y(1) y(2)],'k','LineWidth',2, 'DisplayName', 'Predictions start');
ylabel('error - theta');
xlabel('time');
title('Error between AR and true data');
legend();

