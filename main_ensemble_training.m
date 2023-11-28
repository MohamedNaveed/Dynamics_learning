% ARMA model for dynamics fit. 

clc;clear;
rng(0);

%% define pendulum  
model = pendulum_model(); 


%% ensemble data fit - training
method = 'wDMD'; %'ARMA' | 'wDMD' %ARMA - autoregressive model. wDMD - Window DMD/Hankel DMD.

n_samples = 1;
window = 20; % window / time delayed samples considered for Hankel DMD.
t_span = 10; %total time period
ini_angle = 90; %deg

[A, A_r, error_fit_A, error_fit_Ar, U, S, V] = model_fit(method, model, window, n_samples, t_span, ini_angle);

% A - computed considering all the modes.
% A_r - computed considering 99.999% energy in the modes. 

fprintf('Max training error A = %d \n', max(max(error_fit_A)));
fprintf('Max training error Ar = %d \n', max(max(error_fit_Ar)));
fprintf('Rank of Ar = %d \n', rank(A_r)); 

%% plot singular values of S

figure;
plot(1:size(S,1), diag(S), 'LineWidth',2);
xlabel('index of singular value');
ylabel('singular value');

%% monte_carlo test - testing

t_span = 20;
n_mc_runs = 1;

[error_mean, error_std] = monte_carlo_test(A, method, model,...
                                    window, n_mc_runs, t_span, ini_angle);

%% plot
x_idxs = (window)*model.dt:model.dt:t_span;
y_idxs = window+1:size(error_std,2);
fig = figure;
subplot(2,1,1);
hold on;
title(['Error - ', method, '. window - ', num2str(window)])
plot(x_idxs, error_mean(1,y_idxs),'b','LineWidth',2,'HandleVisibility','off');
%plot(x_idxs, error_mean(1,y_idxs) + error_std(1,y_idxs),'--r','LineWidth',2,...
%    'DisplayName', '1-\sigma');
%plot(x_idxs, error_mean(1, y_idxs) - error_std(1,y_idxs),'--r','LineWidth',2,...
%    'HandleVisibility','off');
%ylim([-1,1]);
ylabel('error - theta');
%legend();

subplot(2,1,2);
hold on;
plot(x_idxs, error_mean(2,y_idxs),'b','LineWidth',2);
%plot(x_idxs, error_mean(2,y_idxs) + error_std(1,y_idxs),'--r','LineWidth',2);
%plot(x_idxs, error_mean(2,y_idxs) - error_std(1,y_idxs),'--r','LineWidth',2);
ylabel('error - theta dot');
xlabel('time');
set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
%saveas(fig,'DMD_w1_tr10_ts10_domain_15.pdf')



%% plot modes
%{
[V_A, D_A] = eig(A);
[V_Ar, D_Ar] = eig(A_r);

%%
eigen_vector_1 = V_A(:,15) + V_A(:,16)
eigen_vector_1 = reshape(eigen_vector_1, [2,10])
figure;
subplot(2,1,1);
plot(1:10,eigen_vector_1(1,:),'LineWidth',2);
subplot(2,1,2);
plot(1:10,eigen_vector_1(2,:),'LineWidth',2);

%%
i = 5;
data_u = reshape(U(:,i),[2,10]);
figure(i);
subplot(2,1,1);
plot(1:10,data_u(1,:),'LineWidth',2);
subplot(2,1,2);
plot(1:10,data_u(2,:),'LineWidth',2);
%}
