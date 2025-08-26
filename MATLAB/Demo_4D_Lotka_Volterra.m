clear, clc, close all;
figpath = './figures/';
addpath('./utils');

%% generate Data
polyorder = 2;
usesine = 0;

n = 4;

x0=[0.5; 0.2; 0.1; 0.5];  % Initial condition

% Integrate
tspan=0:0.01:500;
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(@(t,x) Hyper_LV(t,x),tspan,x0,options);

% figure;
% for i = 1:n
%     subplot(n, 1, i);
%     plot(t, x(:,i));
% end

%% compute Derivative
rng(123456789);
eps = 0.001;
for i=1:length(x)
    dx_free(i,:) = Hyper_LV(0,x(i,:));
end
dx = dx_free + eps*randn(size(dx_free));

% figure;
% for i = 1:n
%     subplot(n, 1, i);
%     plot(t, dx_free(:,i));
%     plot(t, dx(:, i))
% end

%% build library of nonlinear time series 
% Here we focus on pure polynomial basis.
% The variable 'Charset' represents the textual form corresponding to the polynomial basis functions.
[Theta, Charset] = getLibrary(x, polyorder);

%% compute Sparse regression: multiple algorithms
% The sparse optimization problem can also be solved by some other algorithms.
% The only important thing is to obtain the Coefficients Matrix 'Xi'.

% % algorithm 1: sequential thresholded least squares (STLS)
% lambda = 0.025;      % lambda is the threshold value.
% Xi = sparsifyDynamics(Theta,dx,lambda,n);

% % algorithm 2: sequential thresholded ridge regression (STRidge)
% lambda = 0.025;      % lambda is the threshold value.
% Xi = sparsifyDynamics_STRidge(Theta,dx,lambda);

% % algorithm 3: STLS_ABSS
% [Xi, thresholdValue, numIter, history] = sparsifyDynamics_ABSS(Theta, dx);

% algorithm 4: STRidge_ABSS
eta = 1e-6;            % eta is the regularization parameter of ridge regression
[Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, eta);


% % only valid for ABSS kind algorithms
% disp(thresholdValue)
% disp(numIter)

disp(table(Charset, Xi))

%% integrate true and identified systems | The key of FSINDy
func_name = 'test_4D_Hyper_LV';
mg = ModelGenerate(Charset, Xi, func_name); % Automatically convert them to MATLAB function expressions.
% disp(mg); % Readers can display the structure of this variable.

[t_span,x_true]=ode45(@(t,x) Hyper_LV(t,x),tspan,x0,options);   % true model
[~, x_identified] = ode45(func_name,tspan,x0,options);  % approximate

%% evaluate
fpt = zeros(1, n);
epsilon_list = max(abs(x_true)) * 0.1;
for i = 1:n
    epsilon = epsilon_list(i);
    fpt(i) = FPT(t_span, x_true(:, i), x_identified(:, i), epsilon);
    fprintf('X_%d: FPT=%.2f\n',i, fpt(i));
end
dt = mean(diff(t_span));
validPredTime = min(fpt);
validPredInd = floor(validPredTime/dt);
%% Figures
number_label = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

figure('Units', 'centimeters', 'Position', [5 5 17 15]);
t3 = tiledlayout(n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

set(groot, 'defaultAxesFontName', 'Times New Roman',...
    'defaultTextFontName', 'Times New Roman',...
    'defaultAxesLabelFontSizeMultiplier', 1.1,...n
    'defaultAxesTitleFontSizeMultiplier', 1.2);

lineStyle = struct(...
    'True', {'Color', [0.2 0.4 0.8], 'LineWidth', 2.5, 'LineStyle', '-'},... 
    'Identified', {'Color', [0.9 0.4 0.1], 'LineWidth', 2, 'LineStyle', '--'});

for i = 1:n
    ax = nexttile;
    
    plot(t_span, x_true(:,i), lineStyle.True);
    hold on;
    plot(t_span, x_identified(:,i), lineStyle.Identified);
    
    tmp_1 = x_true(:,i);
    tmp_2 = x_identified(:,i);
    fprintf('X_%d: RMSE=%.6f\n', i, error_func(tmp_1, tmp_2, 'rmse'));
    fprintf('X_%d: T=[0, %.2f] RMSE=%.6f\n', i, validPredTime, ...
        error_func(tmp_1(1:validPredInd), tmp_2(1:validPredInd), 'rmse'));

    set(ax, 'LineWidth', 1.5,...          
        'FontSize', 12,...                
        'TickDir', 'out',...              
        'TickLength', [0.015 0.015],...   
        'XMinorTick', 'on',...            
        'YMinorTick', 'on',...
        'Box', 'off');                    
    
    grid(ax, 'on');
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.2;
    ax.GridColor = [0.4 0.4 0.4];
    
    ylabel(['X_',num2str(i)],...
        'FontSize', 12,...
        'FontWeight', 'normal');
    
    title(number_label{i},...
        'FontSize', 12,...
        'HorizontalAlignment', 'left',...
        'Units', 'normalized',...
        'Position', [0.015 0.83],...
        'FontWeight', 'bold');
end

xlabel(t3, 'Time (s)',...
    'FontSize', 12,...
    'FontWeight', 'normal');

leg = legend(ax, {'True','Identified'},...
    'Box', 'off',...
    'Orientation', 'horizontal',...
    'Position', [0.35 0.97 0.5 0.03],...
    'FontSize', 12);
leg.ItemTokenSize = [15,5];