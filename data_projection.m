function [error,mag_error,projections] = data_projection(A,X_data)

[V,D,W] = eig(A);
diag_D = diag(D)

U = inv(V)';

projections = U'*X_data;

c = projections(:,1);
error = zeros(size(projections));
mag_error = zeros(1,size(projections,2));

for k = 2:size(projections,2)
    
    error(:,k) = projections(:,k) - D^(k-1)*c;
    mag_error(k) = norm(error(:,k),2);

end

end