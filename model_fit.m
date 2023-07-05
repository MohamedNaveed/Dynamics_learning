function [A, error_fit, U, S, V] = model_fit(method,model, window, n_samples, t_span)

t_steps = t_span/model.dt;

control = 0;
x = zeros(model.nx, window + 1);

if strcmp(method, 'ARMA')
    
    X = zeros(model.nx*window, n_samples*(t_steps - window + 1));
    Xprime = zeros(model.nx, n_samples*(t_steps - window + 1));

    for n = 1:n_samples

        x0_theta = unifrnd(0,90);
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
    
    X = zeros(model.nx*window, n_samples*(t_steps - window + 1));
    Xprime = zeros(model.nx*window, n_samples*(t_steps - window + 1));

    for n = 1:n_samples

        x0_theta = unifrnd(0,90);
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

    error_fit = Xprime - A*X;
end

end