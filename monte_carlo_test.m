function [error_mean, error_std] = monte_carlo_test(A,method,model, window,...
                    n_mc_runs, t_span, ini_angle)

t_steps = t_span/model.dt;
control = 0;
x = zeros(model.nx, t_steps+1);
error_x1 = zeros(n_mc_runs, t_steps+1);
error_x2 = zeros(n_mc_runs, t_steps+1);

if strcmp(method,'wDMD')
    rng(0);
    x_DMD = zeros(model.nx,t_steps+1);

    for n = 1:n_mc_runs
        x0_theta = unifrnd(0,ini_angle);
        x0 = [deg2rad(x0_theta),0];
        x(:,1) = x0;
        x_DMD(:,1) = x0;
        y_DMD = zeros(model.nx*window,1); %arma state
        y_DMD(model.nx*(window - 1) + 1: model.nx*window) = x0;

        for i = 1:t_steps
    
            x(:,i+1) = pendulum_nl_state_prop(i,x(:,i),control,model); %true

            if i < window
                x_DMD(:,i+1) = x(:,i+1); 
                y_DMD(model.nx*(window - i - 1) + 1:model.nx*(window - i)) = x(:,i+1);
            else
                y_DMD = A*y_DMD;
                x_DMD(:,i+1) = y_DMD(1:model.nx);
            end
        end

        
        x_max = max(x,[],2); % max value for normalization.
        error_x1(n,:) = (x(1,:) - x_DMD(1,:))./x_max(1);
        error_x2(n,:) = (x(2,:) - x_DMD(2,:))./x_max(2);

    end
    
    error_x1_mean = mean(error_x1,1);
    error_x1_std = std(error_x1,0,1);
    error_x2_mean = mean(error_x2,1);
    error_x2_std = std(error_x2,0,1);
    
    error_mean = [error_x1_mean; error_x2_mean];
    error_std = [error_x1_std; error_x2_std];

elseif strcmp(method, 'ARMA')
    rng(0);
    x_arma = zeros(model.nx,t_steps+1);

    for n = 1:n_mc_runs
        x0_theta = unifrnd(0,90);
        x0 = [deg2rad(x0_theta),0];
        x(:,1) = x0;
        x_arma(:,1) = x0;
        y_arma = zeros(model.nx*window,1); %arma state
        y_arma(model.nx*(window - 1) + 1: model.nx*window) = x0;

        for i = 1:t_steps
    
            x(:,i+1) = pendulum_nl_state_prop(i,x(:,i),control,model); %true

            if i < window
                x_arma(:,i+1) = x(:,i+1); 
                y_arma(model.nx*(window - i - 1) + 1:model.nx*(window - i)) = x(:,i+1);
            else
                x_arma(:,i+1) = A*y_arma;
                y_arma = [x_arma(:,i+1);y_arma(1:model.nx*(window-1))];
            end
        end

        x_max = max(x,[],2); % max value for normalization.
        error_x1(n,:) = (x(1,:) - x_arma(1,:))./x_max(1);
        error_x2(n,:) = (x(2,:) - x_arma(2,:))./x_max(2);

    end
    
    error_x1_mean = mean(error_x1,1);
    error_x1_std = std(error_x1,0,1);
    error_x2_mean = mean(error_x2,1);
    error_x2_std = std(error_x2,0,1);
    
    error_mean = [error_x1_mean; error_x2_mean];
    error_std = [error_x1_std; error_x2_std];

end

end