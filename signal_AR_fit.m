% modeling periodic signals using AR model.
clc;clear;

%% define signal
fs = 1; % sampling freq.
dt = 1/fs; % sampling interval
fs_true = 200; % higher frequency to check for aliasing effects. 
dt_true = 1/fs_true;
model.dt = dt;
model.nx = 2;
t_span = 100;
t = 0:dt:t_span;%
%t_true = 0:dt_true:t_span;
%true_model.dt = dt_true;
t_steps = t_span/dt;
f1 = 1;
f2 = 10;
f3 = 11;
f4 = 30;
f5 = 60; %freq out of domain sampling

%x = exp(-0.5*t).*sin(2*pi*(f1).*t);
%x = sin(2*pi*f1*t) + 0.5*sin(2*pi*f2*t);% + 0.5*sin(2*pi*f3*t)+ 0.5*sin(2*pi*f4*t);
%x = sin(2*pi*f1*t)./((t+1).^3);

%x_true = sin(2*pi*f1*t_true) + sin(2*pi*f5*t_true); % to check for
%aliasing effects.
N = 11;

omega_0 = 2*pi/N;
a_actual = [2,1,1,2;1, 0, 0, 1];
a_vec = [1 2;0 1];

x = a_vec(:,1)*(exp(omega_0*1i.*t) + exp(-omega_0*1i.*t))...
    + a_vec(:,2)*(exp(2*omega_0*1i.*t) + exp(-2*omega_0*1i.*t));



x_max = max(x,[],2); 


save_path = "/home/naveed/Documents/Dynamics_learning/plots/signals/test/";

%% plot fft.

fig = fft_signal(x, model,'FFT of True signal');
saveas(fig, save_path + "fft_true.jpg");

%% fit AR model.
window = 2; % window / time delayed samples considered.
n_samples = 20;%1.0/model.dt; % training samples columns of X

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
    
    A_DMD = Xprime*pinv(X); %exact reconstruction of A (includes all the modes)
    
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
    
    A_DMD = Xprime*V(:,1:r)*inv(S(1:r,1:r))*U(:,1:r)'; % A calculated using reduced modes.
    
    A_arma = Xprime_arma*V(:,1:r)*inv(S(1:r,1:r))*U(:,1:r)'; % A calculated using reduced modes.
end



error_training = (Xprime_arma - A_arma*X)./x_max;
RMSE = sqrt(mse(error_training))
%% prediction using A_arma

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
fig = fft_signal(x_arma,model,'FFT of AR predictions');
saveas(fig, save_path + "fft_ar.jpg");
fig = fft_signal(error,model,'FFT of error');
saveas(fig, save_path + "fft_error.jpg");
%%
fig = figure;
plot(t, x, 'k', 'LineWidth',2,'Marker','.','MarkerSize',10, 'DisplayName','Truth');
hold on;
plot(t, x_arma,'color',	[0.8500 0.3250 0.0980], 'LineStyle','-','Marker','.','MarkerSize',1, 'LineWidth',2,'DisplayName','Prediction');
plot(t(1:window+n_samples), x(:,1:window+n_samples),'color','b','LineStyle','-', 'LineWidth',2,'DisplayName','Training data');

legend();
ylabel('signal');
xlabel('time');
title("Window = " + num2str(window) + "; Training samples = " + num2str(n_samples));
saveas(fig, save_path + "pred.jpg");
%%
fig = figure;
hold on;
plot(t, error(1,:),'color',	[0.8500 0.3250 0.0980], 'LineWidth',2, 'HandleVisibility','off'); 
plot(t(1:window+n_samples), error(1,1:window+n_samples), 'color','b','LineStyle','-', 'LineWidth',2,'DisplayName','Training');
y = ylim; % current y-axis limits
x_idx = (window)*model.dt;
plot([x_idx  x_idx],[y(1) y(2)],'k','LineWidth',2, 'DisplayName', 'Predictions start');
ylabel('error');
xlabel('time');
title("Error between AR and true data w=" + num2str(window));
legend();
saveas(fig, save_path + "error.jpg");

%% analyzing A


disp('DMD eigenvalues')
[V,D,W] = eig(A_DMD);
diag_D = diag(D)
frequencies = logm(D)/model.dt;
diag_frequencies = diag(frequencies)

A_DMD_arma = [A_arma; eye(window*model.nx-model.nx) zeros(window*model.nx-model.nx,model.nx)];
disp('AR eigenvalues')
[V,D,W] = eig(A_DMD_arma);
diag_D = diag(D)
frequencies = logm(D)/model.dt;
diag_frequencies = diag(frequencies)

%% analysis using Fourier and Vandermonde matrix

Vand_mat = complex(zeros(window,window)); % vandermonde matrix
freq_idx = [-2,-1,1,2];
time_idx = 2:window+1;%0:window-1;

for row_id = 1:window
    for col_id = 1:window

        Vand_mat(row_id,col_id) = exp(1i*omega_0*freq_idx(row_id)*time_idx(col_id));

    end
end

x_pred_modes = complex(zeros(window,1));

pred_idx = window+2;
for k = 1:window
    x_pred_modes(k) = exp(1i*omega_0*freq_idx(k)*pred_idx);
end

%% check if prediction are correct
x_pred_from_a = a_actual*Vand_mat;

arma_par_fourier = Vand_mat\x_pred_modes;

est_a = x(:,pred_idx-window+1:pred_idx)/Vand_mat;

%% checking if the scalar arma from fourier is equivalent to vector arma
A_arma_fourier = zeros(model.nx,model.nx*window);

for i=1:window
    A_arma_fourier(:,(i-1)*model.nx+1:i*model.nx) = arma_par_fourier(window-i+1)*eye(model.nx);
end

x_arma_fourier = zeros(model.nx,t_steps+1);

y_arma_fourier = zeros(model.nx*window,1); %observable

%initial condition
for w = 1:window
    y_arma_fourier(model.nx*(window - w) + 1:model.nx*(window - w + 1)) = x(:,w);
end

x_arma_fourier(:,1:window) = x(:,1:window); %assume same for the first few steps. 

for i = window:t_steps

    x_arma_fourier(:,i+1) = A_arma_fourier*y_arma_fourier; 
    y_arma_fourier = [x_arma_fourier(:,i+1);y_arma_fourier(1:model.nx*(window-1))];
end

%% error 

error_fourier = (x - x_arma_fourier)./x_max;



