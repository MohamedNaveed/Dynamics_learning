function [A, A_r, error_fit, error_fit_Ar, U, S, V] = model_fit(method,model, window, n_samples, t_span, ini_angle)

t_steps = t_span/model.dt;

control = 0;
x = zeros(model.nx, window + 1);

if strcmp(method, 'ARMA')
    rng(0);
    X = zeros(model.nx*window, n_samples*(t_steps - window + 1));
    Xprime = zeros(model.nx, n_samples*(t_steps - window + 1));

    for n = 1:n_samples

        if n_samples > 1
            x0_theta = unifrnd(0,ini_angle); %if randomized samples
        else
            x0_theta = ini_angle;
        end
        x0 = [deg2rad(x0_theta),0];
        x(:,1) = x0;
        y_arma = zeros(model.nx*window,1); %arma state
        y_arma(model.nx*(window - 1) + 1: model.nx*window) = x0;

        for i = 1:t_steps

            x(:,i+1) = model.state_propagate(i,x(:,i),control,model); %true
            
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
    fprintf('Rank of X matrix: %d \n', rank(X));
    fprintf('Number of rows in X: %d \n', size(X,1));

    A = Xprime*pinv(X);

    error_fit = Xprime - A*X;

elseif strcmp(method, 'wDMD')
    rng(0);
    X = zeros(model.nx*window, n_samples*(t_steps - window + 1));
    Xprime = zeros(model.nx*window, n_samples*(t_steps - window + 1));

    for n = 1:n_samples
        
        if n_samples > 1
            x0_theta = unifrnd(0,ini_angle); %if randomized samples
        else
            x0_theta = ini_angle;
        end
        
        x0 = [deg2rad(x0_theta),0];
        x(:,1) = x0;
        y = zeros(model.nx*window,1); %concatenated state / stacked
        y(model.nx*(window - 1) + 1: model.nx*window) = x0;
        
        % building data matrix. 
        for i = 1:t_steps

            x(:,i+1) = model.state_propagate(i,x(:,i),control,model); %true (simulating dynamics)
            
            if i < window
            
                y(model.nx*(window - i - 1) + 1:model.nx*(window - i)) = x(:,i+1);
            
            else
                idx = (n-1)*(t_steps - window + 1) + (i - window) + 1;
                X(:, idx) = y;  % build data matrix
                 
                y = [x(:,i+1);y(1:model.nx*(window-1))];
                Xprime(:,idx) = y;

            end
            
        end
        
        
        
    end
    fprintf('Rank of X matrix: %d \n', rank(X));
    fprintf('Number of rows in X: %d \n', size(X,1));

    [U, S, V] = svd(X); %SVD of the data matrix 
    
    A = Xprime*pinv(X); %exact reconstruction of A (includes all the modes) 
    
    thresh = 0.99999; % threshold for selecting the modes 
    diag_S = diag(S);
    sum_S = sum(diag_S); % sum of all the singular values. 
    
    % finding the number of modes to meet the threshold (thresh) value.
    for S_i = 1:length(diag_S)
       
        sum_S_i = sum(diag_S(1:S_i)); % partial sum of all the singular values. 
        
        if sum_S_i/sum_S >= thresh
            r = S_i;
            break;
        end      
    end
    
    A_r = Xprime*V(:,1:r)*inv(S(1:r,1:r))*U(:,1:r)'; % A calculated using reduced modes.
    
    error_fit_Ar = Xprime - A_r*X; % training error using Ar
    error_fit = Xprime - A*X;  %training error using A
end

end