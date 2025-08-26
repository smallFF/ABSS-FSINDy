function [Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, DX, lambda, MRE_threshold) 
if nargin<4
    MRE_threshold = 0.01;
end
objective = @(A, x, b) 1/2*norm(A*x-b).^2;
A = Theta;
[~, p] = size(A);
dim = size(DX, 2);

AtA = Theta'*Theta;
I = eye(size(AtA));
Xi = (AtA + lambda*I) \ (A'*DX);


numIter = zeros(1, dim);
thresholdValue = zeros(1, dim);

history.objval = cell(1, dim);
history.MRE = cell(1, dim);
history.FValue = cell(1, dim);
history.deletedValue = cell(1, dim); % Record redundant coefficients.
history.Ftest = cell(1, dim);
history.indmat = cell(1, dim);
for i = 1:dim
    b = DX(:, i);
    x_init = Xi(:, i);

    x = x_init;
    
    R = [ ];
    ind = logical(1:p);
    history.indmat{i} = zeros(p, p);
    history.indmat{i}(:, 1) = ind;
    for k = 1:p-1
        ind_old = ind;

        x_old = x;
        x_min = min(abs(x(ind)));
        threshold = x_min + eps;

        ind = abs(x)>threshold;
        x(~ind) = 0;

        A_tmp = A(:, ind);
        AtA_tmp = A_tmp'*A_tmp;
        I_tmp = eye(size(AtA_tmp));
        x(ind) = (AtA_tmp + lambda*I_tmp) \ (A_tmp'*b);

        history.indmat{i}(:, 1+k) = ind;
        
        % redord history information
        history.objval{i}(k) = objective(A, x, b);
        obj_old = objective(A, x_old, b);
        obj_new = objective(A, x, b);
        history.MRE{i}(k) = calModifiedRelativeError(obj_old, obj_new);

        history.FValue{i}(k) = ((obj_new - obj_old)/obj_old)*k;
       
        R = [R, threshold];
        history.deletedValue{i}(k) = threshold;

        if history.MRE{i}(k) > MRE_threshold
            numIter(i) = k-1;
            thresholdValue(i) = max( R(1:k-1) );
            x = x_old;
            break;
        end
        
        % base model (only one parameter)
        if sum(ind) == 1
            numIter(i) = k;
            thresholdValue(i) = max( R(1:k) );
        end
    end
    
    Xi(:, i) = x;
end
