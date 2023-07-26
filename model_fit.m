function [A, A_r, error_fit, error_fit_Ar, U, S, V] = model_fit(method,model, window, n_samples, t_span, ini_angle)

t_steps = t_span/model.dt;

control = 0;
x = zeros(model.nx, window + 1);

if strcmp(method, 'ARMA')
    rng(0);
    X = zeros(model.nx*window, n_samples*(t_steps - window + 1));
    Xprime = zeros(model.nx, n_samples*(t_steps - window + 1));

    for n = 1:n_samples

        x0_theta = unifrnd(0,ini_angle);
        x0 = [deg2rad(x0_theta),0];
        x(:,1) = x0;
        y_arma = zeros(model.nx*window,1); %arma state
        y_arma(model.nx*(window - 1) + 1: model.nx*window) = x0;

        for i = 1:t_steps

            x(:,i+1) = pendulum_nl_state_prop(i,x(:,i),control,model); %true
            
            if i < window
            
                y_arma(model.nx*(window - i - 1) + 1:model.nx*(window - i)) = x(:,i+1);
            
            else
                idx = (n-1)*(t_steps - window + 1) + (i - window) + 1;
                X(:, idx) = y_arma;  % build data matrix
                Xprime(:,idx) = x(:,i+1); 
                y_arma = [x(:,i+1);y_arma(1:model.nx*(window-1))];

            end
            
        end
        
        
        
    end
    disp('rank of X =')
    rank(X)
    disp('number of rows of X');
    size(X,1)

    A = Xprime*pinv(X);

    error_fit = Xprime - A*X;

elseif strcmp(method, 'wDMD')
    rng(0);
    X = zeros(model.nx*window, n_samples*(t_steps - window + 1));
    Xprime = zeros(model.nx*window, n_samples*(t_steps - window + 1));

    for n = 1:n_samples

        x0_theta = unifrnd(0,ini_angle);
        x0 = [deg2rad(x0_theta),0];
        x(:,1) = x0;
        y_arma = zeros(model.nx*window,1); %arma state
        y_arma(model.nx*(window - 1) + 1: model.nx*window) = x0;

        for i = 1:t_steps

            x(:,i+1) = pendulum_nl_state_prop(i,x(:,i),control,model); %true
            
            if i < window
            
                y_arma(model.nx*(window - i - 1) + 1:model.nx*(window - i)) = x(:,i+1);
            
            else
                idx = (n-1)*(t_steps - window + 1) + (i - window) + 1;
                X(:, idx) = y_arma;  % build data matrix
                 
                y_arma = [x(:,i+1);y_arma(1:model.nx*(window-1))];
                Xprime(:,idx) = y_arma;

            end
            
        end
        
        
        
    end
    disp('rank of X =')
    rank(X)
    disp('number of rows of X');
    size(X,1)
    [U, S, V] = svd(X);
    A = Xprime*pinv(X);
    thresh = 0.99999;
    diag_S = diag(S);
    sum_S = sum(diag_S);
    
    for S_i = 1:length(diag_S)
       
        sum_S_i = sum(diag_S(1:S_i));
        
        if sum_S_i/sum_S >= thresh
            r = S_i;
            break;
        end      
    end
    
    A_r = Xprime*V(:,1:r)*inv(S(1:r,1:r))*U(:,1:r)';
    
    error_fit_Ar = Xprime - A_r*X;
    error_fit = Xprime - A*X;
end

end