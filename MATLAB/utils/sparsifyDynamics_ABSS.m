function [Xi, thresholdValue, numIter, history] = sparsifyDynamics_ABSS(Theta, DX, MRE_threshold)
if nargin<3
    MRE_threshold = 0.05;
end
objective = @(A, x, b) 1/2*norm(A*x-b).^2;
A = Theta;
[~, p] = size(A);
dim = size(DX, 2);
Xi = zeros(p, dim);

% D = diag(vecnorm(A));
% A = A / D;
% diag_L2_norm = L2_norm(A);
% A = A / diag_L2_norm;

numIter = zeros(1, dim);
thresholdValue = zeros(1, dim);

history.objval = cell(1, dim);
history.MRE = cell(1, dim);
history.FValue = cell(1, dim);
history.deletedValue = cell(1, dim); % Record redundant coefficients.
% history.AdjRS = cell(1, dim);
history.indmat = cell(1, dim);
for i = 1:dim
    b = DX(:, i);
    x_init = A \ b;

    x = x_init;
    
    R = [ ];
    ind = logical(1:p);
    % fprintf("%d | initial index:\n", i);
    % disp(ind);
    history.indmat{i} = zeros(p, p);
    history.indmat{i}(:, 1) = ind;
    for k = 1:p-1
        x_old = x;
        x_min = min(abs(x(ind)));
        threshold = x_min + eps;

        ind = abs(x)>threshold;
        x(~ind) = 0;
        x(ind) = A(:,ind) \ b;

        % fprintf("iter=%d\n", k);
        % disp(ind);
        history.indmat{i}(:, 1+k) = ind;
        
        % redord history information
        history.objval{i}(k) = objective(A, x, b);
        obj_old = objective(A, x_old, b);
        obj_new = objective(A, x, b);
        history.MRE{i}(k) = calModifiedRelativeError(obj_old, obj_new);
        history.FValue{i}(k) = ((obj_new - obj_old)/obj_old)*k;
       
        R = [R, threshold];
        history.deletedValue{i}(k) = threshold;

        % ms = ModelSelection(A, x_old, b);
        % % RSS;    % RSS: Residual Sum of Squares
        % % MSE;    % MSE: Mean Square Error
        % % TSS;    % TSS: Total Sum of Squares
        % % RSE;    % RSE: Redidual Standard Error
        % % RS;     % RS: R^2
        % % Cp;
        % % AIC;    % AIC: Akaike Information Criterion
        % % BIC;    % BIC: Bayesian Information Criterion
        % % AdjRS;  % AdjRS: Adjusted R^2
        % history.RSS(k, i) = ms.RSS;
        % history.MSE(k, i) = ms.MSE;
        % history.TSS(k, i) = ms.TSS;
        % history.RSE(k, i) = ms.RSE;
        % history.RS(k, i) = ms.RS;
        % history.Cp(k, i) = ms.Cp;
        % history.AIC(k, i) = ms.AIC;
        % history.BIC(k, i) = ms.BIC;
        % history.RS{i}(k) = ms.RS;
        % history.AICc(k, i) = ms.AICc;

%         disp(eval("ms.Cp"));
%         if history.MRE(k, i) > 0.05 || sum(ind)<=1
        % fprintf("%d : history.MRE=%s, ind=%s\n", i, num2str(history.MRE(k,i)), num2str(sum(ind)))
%         if sum(ind)<=1

        % abrupt change
        if history.MRE{i}(k) > MRE_threshold
        % if history.FValue{i}(k) > 1
        % if history.MRE{i}(k) > 0.01
        % if sum(ind)>1
            numIter(i) = k-1;

            % thresholdValue(i) = max(history.deletedValue{i}(1:k-1));
            thresholdValue(i) = max( R(1:k-1) );
            x = x_old;

            % ms = ModelSelection(A, x_old, b);
            % % RSS;    % RSS: Residual Sum of Squares
            % % MSE;    % MSE: Mean Square Error
            % % TSS;    % TSS: Total Sum of Squares
            % % RSE;    % RSE: Redidual Standard Error
            % % RS;     % RS: R^2
            % % Cp;
            % % AIC;    % AIC: Akaike Information Criterion
            % % BIC;    % BIC: Bayesian Information Criterion
            % % AdjRS;  % AdjRS: Adjusted R^2
            % history.RSS(k, i) = ms.RSS;
            % history.MSE(k, i) = ms.MSE;
            % history.TSS(k, i) = ms.TSS;
            % history.RSE(k, i) = ms.RSE;
            % history.RS(k, i) = ms.RS;
            % history.Cp(k, i) = ms.Cp;
            % history.AIC(k, i) = ms.AIC;
            % history.BIC(k, i) = ms.BIC;
            % history.RS{i}(k) = ms.RS;
            % 
            % history.AICc(k, i) = ms.AICc;
            break;
        end
        
        % base model (only one parameter)
        if sum(ind) == 1
        % if sum(ind)>1
            numIter(i) = k;
            % thresholdValue(i) = max(history.deletedValue{i}(1:k));
            thresholdValue(i) = max( R(1:k) );
            % x
            % x = x_old;

            % ms = ModelSelection(A, x_old, b);
            % % RSS;    % RSS: Residual Sum of Squares
            % % MSE;    % MSE: Mean Square Error
            % % TSS;    % TSS: Total Sum of Squares
            % % RSE;    % RSE: Redidual Standard Error
            % % RS;     % RS: R^2
            % % Cp;
            % % AIC;    % AIC: Akaike Information Criterion
            % % BIC;    % BIC: Bayesian Information Criterion
            % % AdjRS;  % AdjRS: Adjusted R^2
            % history.RSS(k, i) = ms.RSS;
            % history.MSE(k, i) = ms.MSE;
            % history.TSS(k, i) = ms.TSS;
            % history.RSE(k, i) = ms.RSE;
            % history.RS(k, i) = ms.RS;
            % history.Cp(k, i) = ms.Cp;
            % history.AIC(k, i) = ms.AIC;
            % history.BIC(k, i) = ms.BIC;
            % history.AdjRS(k, i) = ms.AdjRS;
            % 
            % history.AICc(k, i) = ms.AICc;
            % break;
        end
    end
    
    Xi(:, i) = x;
end

% disp(history)
% Xi = diag_L2_norm \ Xi;
% Xi = D \ Xi;


% mean_objval = history.objval./mean(history.objval);
% 
% negligible_error = 5e-3;
% figure;
% for i = 1:dim
%     subplot(dim, 1, i);
%     plot(log10(mean_objval(:, i)), '.-', 'LineWidth',1.5);
%     xlabel('Iterations', 'FontSize',13)
%     ylabel('Objective Value', 'FontSize',13)
% 
%     criterion_value = log(mean_objval(:, i));
% 
%     plot(criterion_value);
%     xlabel('Iterations', 'FontSize',13)
%     ylabel('Objective Value', 'FontSize',13)
%     criterion_min = min(criterion_value);
%     criterion_max = max(criterion_value);
%         % fprintf('Fig(%d), %f\n', 100+i, criterion_max-criterion_min);
%     ylim([criterion_min-negligible_error, criterion_max+negligible_error]);
% end

