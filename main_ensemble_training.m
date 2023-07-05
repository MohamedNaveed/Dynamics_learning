% ARMA model for dynamics fit. 

clc;clear;
rng(0);
%% define pendulum  
model = pendulum_model(); 


%% ensemble data fit - training
method = 'wDMD'; %'ARMA' | 'wDMD'
n_samples = 100;
window = 10;
t_span = 6; %total time period

[A, error_fit, U, S, V] = model_fit(method, model, window, n_samples, t_span);

fprintf('Max training error = %d \n', max(max(error_fit)));

%% plot singular values of S

figure;
plot(1:size(S,1), diag(S), 'LineWidth',2);

%% monte_carlo test - testing

t_span = 10;
n_mc_runs = 100;

[error_mean, error_std] = monte_carlo_test(A, method, model,...
                                    window, n_mc_runs, t_span);

%% plot
x_idxs = (window)*model.dt:model.dt:t_span;
y_idxs = window+1:size(error_std,2);
fig = figure;
subplot(2,1,1);
hold on;
title([method, '. window - ', num2str(window)])
plot(x_idxs, error_mean(1,y_idxs),'b','LineWidth',2,'HandleVisibility','off');
plot(x_idxs, error_mean(1,y_idxs) + error_std(1,y_idxs),'--r','LineWidth',2,...
    'DisplayName', '1-\sigma');
plot(x_idxs, error_mean(1, y_idxs) - error_std(1,y_idxs),'--r','LineWidth',2,...
    'HandleVisibility','off');

ylabel('error - theta');
legend();
subplot(2,1,2);
hold on;
plot(x_idxs, error_mean(2,y_idxs),'b','LineWidth',2);
plot(x_idxs, error_mean(2,y_idxs) + error_std(1,y_idxs),'--r','LineWidth',2);
plot(x_idxs, error_mean(2,y_idxs) - error_std(1,y_idxs),'--r','LineWidth',2);

ylabel('error - theta dot');
xlabel('time');

saveas(fig,'DMD_w10_tr6_ts10.jpg')
