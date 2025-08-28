clear, clc, close all;
figpath = './figures/';
addpath('./utils');
isOutputToFile = 0; % Switch for outputting figures in PDF format.

%% generate Data
polyorder = 2;
usesine = 0;

F = 8;

% n = 4, 10, 15
n = 15;

x0 = 8*ones(1, n);
x0(1) = 1;

% Integrate
tspan=[0:0.01:15];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));

func_name_true = generate_lorenz96(n, F);

[t,x]=ode45(func_name_true,tspan,x0,options);

%%
figure('Units', 'centimeters', 'Position', [2 2 17 22.5]);
for i = 1:n
    subplot(n, 1, i);
    plot(t, x(:, i));
    ylabel(['X_{',num2str(i),'}']);
end

%% compute Derivative
rng(123456789);
dx = zeros(N, n);
eps = 1;
for i=1:length(x)
    dx(i,:) = eval([func_name_true,'(0,x(i,:))']);
end
dx = dx + eps*randn(size(dx));

%% build library of nonlinear time series
[Theta, Charset] = getLibrary(x, polyorder);

%% compute Sparse regression: STRidge_ABSS
% [Xi, thresholdValue, numIter, history] = sparsifyDynamics_ABSS(Theta, dx);
[Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);

% disp(thresholdValue)
% disp(numIter)

disp(table(Charset, Xi))
% For these results that are too long to display, we can output them to a file.
% For example
data_file = 'Charset_Xi.csv';
writetable(table(Charset, Xi), data_file);
% delete(data_file); % uncomment this line and run this single line to quick delete the data file.

fprintf("%10s | %-9s | %-10s\n","State", 'Iteration', 'Threshold')
for i = 1:n
    fprintf("%10s | %-9s | %-.6f\n", ['DX_', num2str(i)], num2str(numIter(i)), thresholdValue(i));
end

baseFileName = mfilename;

figure('Units', 'centimeters', 'Position', [2 2 17 22.5]);
t1 = tiledlayout(n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
fontsize = 13;
for i = 1:n
    eval(['ax',num2str(i),'=nexttile;']);
    yyaxis left;
    y1 = history.objval{i}(:);
    plot(y1, '-', 'LineWidth',1.5, 'Marker', '*', 'MarkerSize', 10);
    ylim( [yrange_extend(y1)] );
    
    yyaxis right;
    y2 = history.MRE{i}(:);
    plot(y2, '.-', 'LineWidth',1.5, 'Marker', 'o', 'MarkerSize', 10);
    ylim( [yrange_extend(y2)] );

    grid minor;
    set(gca, 'LineWidth',1.2, 'FontSize', 12, 'FontName', 'times');
    ylabel(['X_{',num2str(i),'}'], 'FontSize',12);

end
xlabel('Iterations', 'FontSize',12)
legend(ax1, 'Objective Value', 'Modified Relative Error', 'Location', 'northoutside', 'Orientation', 'horizontal');
legend(ax1, 'boxoff');

coefTrue = cell(1, n);
for i = 1:n
    coefTrue{i} = [1 -1 -1 8];
end

xlim_range = [1e-5 1e2];
figure;
Xi_ls = Theta \ dx;
for i = 1:n
    x1 = abs(Xi_ls(:, i));
    y1 = -i*ones(size(x1));
    x2 = abs(coefTrue{i});
    y2 = -i*ones(size(x2));

    semilogx(x1, y1, 'r.', 'MarkerSize', 10);
    xlim(xlim_range);

    hold on;

    semilogx(x2, y2, 'b^', 'MarkerSize', 10);

    set(gca, 'LineWidth',1.2, 'FontSize', fontsize, 'FontName', 'Times');

    text(3.5*xlim_range(1)/10, -i, ['X_{',num2str(i),'}'], 'FontSize', fontsize);
    set(gca, 'yTick', [], 'FontName', 'Times');

    if i == 1
        grid minor;
    end
end
xlabel('Absolute value of coefficients (logarithmic coordinate)', 'FontSize',fontsize)
legend(gca, 'Initial Least-square Estimation', 'Key Features', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 12);
legend(gca, 'boxoff');
ylim([-(n+1) 0]);
if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_n_',num2str(n),'_a.pdf'],'ContentType','vector', 'BackgroundColor','none');
end

figure;
P = zeros(n, 2);
for i = 1:n
    x1 = history.deletedValue{i}(1:numIter(i));
    y1 = -i*ones(size(x1));
    x2 = abs(coefTrue{i});
    y2 = -i*ones(size(x2));

    lowerBoundary = max(x1);
    upperBoundary = min(x2);

    P(i, 1) = lowerBoundary;
    P(i, 2) = upperBoundary;

    if lowerBoundary < upperBoundary
        str_Pi = sprintf("$P_{%d}=(%.6f, %.6f)$", i, lowerBoundary, upperBoundary);
    else
        str_Pi = "$\P_{i}=\emptyset$";
    end

    semilogx(x1, y1, 'ro');
    xlim(xlim_range);
    
    hold on;

    line([lowerBoundary, upperBoundary], [-i, -i], 'Color', 'black', 'LineWidth', 2);

    semilogx(x2, y2, 'b^');
    xlim(xlim_range);  
    
    set(gca, 'LineWidth',1.2, 'FontSize', fontsize, 'FontName', 'Times');

    text(3.5*xlim_range(1)/10, -i, ['X_{',num2str(i),'}'], 'FontSize', fontsize);

    text(lowerBoundary, -i+0.25, str_Pi, 'FontName', 'Times', 'FontWeight', 'Bold', 'FontSize', 10, 'Interpreter','latex');

    set(gca, 'yTick', [], 'FontName', 'Times');

    if i == 1
        hold on;
        grid minor;
    end
end
xlabel('Absolute value of coefficients (logarithmic coordinate)', 'FontSize', fontsize)
legend(gca, 'Redundant Features', 'Threshold Interval', 'Key Features', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 12);
legend(gca, 'boxoff');
ylim([-(n+1) 0]);

Pi_min = max(P(:, 1));
Pi_max = min(P(:, 2));
if Pi_min< Pi_max
    str_Pi = sprintf("$\\cap P_{i}=(%.6f, %.6f)$", Pi_min, Pi_max);
else
    str_Pi = "$\cap P_{i}=\emptyset$";
end

line([Pi_min, Pi_min], [-(n+1), 0], 'lineStyle', '--', 'HandleVisibility','off');
line([Pi_max, Pi_max], [-(n+1), 0], 'lineStyle', '--', 'HandleVisibility','off');

text(Pi_min, -n-(0.25+0.21*n/15), str_Pi, 'FontName', 'Times', 'FontWeight', 'Bold', 'FontSize', 10, 'Interpreter','latex');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_n_',num2str(n),'_b.pdf'],'ContentType','vector', 'BackgroundColor','none');
end

%% integrate true and identified systems

func_name = ['test_lorenz96_n_',num2str(n)];
mg = ModelGenerate(Charset, Xi, func_name);

[t_span,x_true]=ode45(func_name_true,tspan,x0,options);   % true model
[~, x_identified] = ode45(func_name,tspan,x0,options);  % approximate
delete([func_name,'.m']);
delete([func_name_true,'.m']);

fpt = zeros(1, n);
epsilon = max( max(abs(x_true)) ) * 0.05;
for i = 1:n
    fpt(i) = FPT(t_span, x_true(:, i), x_identified(:, i), epsilon);
    fprintf('X_%d: FPT=%.2f\n',i, fpt(i));
end
dt = mean(diff(t_span));
validPredTime = min(fpt);
validPredInd = floor(validPredTime/dt);
%% Figures
figure('Position',[100, 100, 1000, 600]);
t3 = tiledlayout(5, 3);
t3.TileSpacing = 'compact';
t3.Padding = 'compact';

len = floor(length(t_span)/3);

for i = 1:n
    eval(['ax',num2str(i),'=nexttile;']);
    tmp_1 = x_true(:, i);
    plot(t_span(1:len), tmp_1(1:len), 'r-', 'LineWidth',2);
    hold on;
    ylim( [yrange_extend(tmp_1(1:len))] );

    tmp_2 = x_identified(:,i);
    plot(t_span(1:len), tmp_2(1:len), 'b--', 'LineWidth',2);
    ylim( [yrange_extend(tmp_2(1:len))] );

    grid on;
    set(gca, 'LineWidth',1.5, 'FontSize', 14, 'FontName', 'Times');
    ylabel(['X_{',num2str(i),'}'], 'FontSize',14);

    fprintf('X_%d: RMSE=%.6f\n', i, error_func(tmp_1, tmp_2, 'rmse'));
    fprintf('X_%d: T=[0, %.2f] RMSE=%.6f\n', i, validPredTime, ...
        error_func(tmp_1(1:validPredInd), tmp_2(1:validPredInd), 'rmse'));
end
xlabel(t3, 'Time', 'FontSize',14);
leg = legend(ax1, 'True', 'Identified', 'Location', 'northoutside', 'Orientation', 'horizontal');
legend(ax1, 'boxoff');
leg.Layout.Tile = 'north';
set(leg, 'FontSize', 15);

if isOutputToFile
    exportgraphics(t3,[figpath, baseFileName,'_n_',num2str(n),'_c.pdf'],'ContentType','vector');
end