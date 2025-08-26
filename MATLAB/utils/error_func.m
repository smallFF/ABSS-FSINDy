function r = error_func(s1, s2, name)
if strcmp(name, 'mse') % mean square error
    func = @(a, b) mean((a-b).^2);
elseif strcmp(name, 'rmse') % root mean square error
    func = @(a, b) sqrt(mean((a-b).^2));
elseif strcmp(name, 'mae') % mean absolute error
    func = @(a, b) mean(sum(abs(a-b)));
elseif strcmp(name, 'rmsre') % root mean squared relative error | sqrt( mean( ((a-b)./a).^2) )
    func = @(a, b) sqrt(mean( ( (a-b)./(1 + abs(a))).^2 ));
elseif strcmp(name, 'rrmse') % relative root mean square error
    func = @(a, b) sqrt(mean((a-b).^2)) ./ mean(a) * 100;
    % disp(sqrt(mean((s1-s2).^2)));
    % disp(mean(s1));
elseif strcmp(name, 'r2') % R square
    func = @(a, b) 1 - (sum( (a-b).^2 ) / sum( (mean(a) - a).^2 ));
end
% fprintf('name=%s\n', name);
r = func(s1, s2);