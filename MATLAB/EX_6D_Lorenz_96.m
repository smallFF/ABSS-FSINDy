clear, clc, close all;
figpath = './figures/';
addpath('./utils');
isOutputToFile = 0; % Switch for outputting figures in PDF format.

%% generate Data
% polyorder = 5;
polyorder = 2;
usesine = 0;

n = 6;
F = 8;

func_name_true = generate_lorenz96(n, F);

x0 = [1 8 8 8 8 8];

% Integrate
tspan=[0:0.01:15];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(func_name_true,tspan,x0,options);

%%
% figure('Units', 'centimeters', 'Position', [2 2 17 22.5]);
% for i = 1:n
%     subplot(n, 1, i);
%     plot(t, x(:, i));
%     ylabel(['X_',num2str(i)]);
% end

%% compute Derivative
rng(123456789);
eps = 1;
% eps = 0;
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

baseFileName = mfilename;

number_label = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

figure('Units', 'centimeters', 'Position', [2 2 17 22.5]);
t1 = tiledlayout(n, 1);
t1.TileSpacing = 'compact';
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

    disp(yrange_extend(y2));
    grid minor;
    set(gca, 'LineWidth',1.2, 'FontSize', 12, 'FontName', 'Arial');
    ylabel(['X_',num2str(i)], 'FontSize',12);
    
    title(number_label{i}, 'Units', 'normalized', 'Position', [0.04, 0.7, 0], 'FontWeight', 'bold');
end
xlabel('Iterations', 'FontSize',12)
legend(ax1, 'Objective Value', 'Modified Relative Error', 'Location', 'northoutside', 'Orientation', 'horizontal');
legend(ax1, 'boxoff');

if isOutputToFile
    exportgraphics(t1,[figpath, baseFileName,'_MRE.pdf'],'ContentType','vector');
end

%%
% figure;
% for i = 1:n
%     subplot(n, 1, i);
%     plot(history.FValue{i}(:), '.-', 'LineWidth',1.5, 'Marker', 'o', 'MarkerSize', 10);
% end

coefTrue = cell(1, n);
for i = 1:n
    coefTrue{i} = [1 -1 -1 8];
end

xlim_range = [1e-6 1e2];
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
    xlim(xlim_range);

    text(4*xlim_range(1)/10, -i, ['X_{',num2str(i),'}']);
    set(gca, 'yTick', []);

    if i == 1
        grid minor;
    end
end
xlabel('Absolute value of coefficients (logarithmic coordinate)', 'FontSize',12)
legend(gca, 'Initial Least-square Estimation', 'Key Features', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 11);
legend(gca, 'boxoff');
ylim([-(n+1) 0]);
if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_a.pdf'],'ContentType','vector', 'BackgroundColor','none');
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
    
    
    set(gca, 'LineWidth',1.2, 'FontSize', 12, 'FontName', 'Arial');

    text(4*xlim_range(1)/10, -i, ['X_{',num2str(i),'}']);

    text(lowerBoundary, -i+0.25, str_Pi, 'FontName', 'Times', 'FontWeight', 'Bold', 'FontSize', 12, 'Interpreter','latex');

    set(gca, 'yTick', []);

    if i == 1
        hold on;
        grid minor;
    end
end
xlabel('Absolute value of coefficients (logarithmic coordinate)', 'FontSize',12)
legend(gca, 'Redundant Features', 'Threshold Interval', 'Key Features', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 11);
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

text(Pi_min, -n-0.25, str_Pi, 'FontName', 'Times', 'FontWeight', 'Bold', 'FontSize', 12, 'Interpreter','latex');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_b.pdf'],'ContentType','vector', 'BackgroundColor','none');
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
number_label = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

figure('Units', 'centimeters', 'Position', [2 2 17 22.5]);
t3 = tiledlayout(n, 1);
t3.TileSpacing = 'compact';
for i = 1:n
    eval(['ax',num2str(i),'=nexttile;']);
    tmp_1 = x_true(:, i);
    plot(t_span, tmp_1, 'r-', 'LineWidth',1.5);
    hold on;
    ylim( [yrange_extend(tmp_1)] );

    tmp_2 = x_identified(:,i);
    plot(t_span, tmp_2, 'b--', 'LineWidth',1.5);
    ylim( [yrange_extend(tmp_2)] );

    grid minor;
    set(gca, 'LineWidth',1.2, 'FontSize', 12, 'FontName', 'Arial');
    ylabel(['X_',num2str(i)], 'FontSize',12);
    
    title(number_label{i}, 'Units', 'normalized', 'Position', [0.04, 0.7, 0], 'FontWeight', 'bold');

    fprintf('X_%d: RMSE=%.6f\n', i, error_func(tmp_1, tmp_2, 'rmse'));
    fprintf('X_%d: T=[0, %.2f] RMSE=%.6f\n', i, validPredTime, ...
        error_func(tmp_1(1:validPredInd), tmp_2(1:validPredInd), 'rmse'));
end
xlabel('Time', 'FontSize',12)
legend(ax1, 'True', 'Identified', 'Location', 'northoutside', 'Orientation', 'horizontal');
legend(ax1, 'boxoff');

if isOutputToFile
    exportgraphics(t3,[figpath, baseFileName,'_c.pdf'],'ContentType','vector');
end