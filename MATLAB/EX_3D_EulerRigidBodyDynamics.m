clear, clc, close all;
figpath = './figures/';
addpath('./utils');
isOutputToFile = 0; % Switch for outputting figures in PDF format.
baseFileName = mfilename;

%% generate Data
polyorder = 3;

n = 3;

x0=[1; 1; 1];  % Initial condition

I1 = 1;
I2 = 2;
I3 = 3;

erbd = @(t, x) [
    -(I3-I2)/I1*x(2)*x(3);
    -(I1-I3)/I2*x(3)*x(1);
    -(I2-I1)/I3*x(1)*x(2);
];


% Integrate
tspan=[0:0.01:50];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(@(t,x) erbd(t,x),tspan,x0,options);

% figure;
% for i = 1:n
%     subplot(n, 1, i);
%     plot(t, x(:, i));
% end

%% compute Derivative
rng(123456789);
eps = 0.01;
for i=1:length(x)
    dx(i,:) = erbd(0,x(i,:));
end
dx = dx + eps*randn(size(dx));

%% pool Data  (i.e., build library of nonlinear time series)
[Theta, Charset] = getLibrary(x, polyorder);

% qualitative analysis
correlation_analysis(Theta, Charset)

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_qualitative_analysis.pdf'],'ContentType','vector', 'BackgroundColor', 'none');
end

% quantitative analysis
figure('Units', 'centimeters', 'Position', [3 3 17 15]);
[U, S, V] = svd(Theta);
semilogy(diag(S), 'o-', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Number', 'FontSize',15);
ylabel('Singular Value [log]', 'FontSize',15);
grid minor
set(gca, 'LineWidth', 1.5, 'FontName', 'Times', 'FontSize', 15);

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_quantitative_analysis.pdf'],'ContentType','vector', 'BackgroundColor', 'none');
end
%% compute Sparse regression: STRidge_ABSS
[Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);

% disp(thresholdValue)
% disp(numIter)
disp(table(Charset, Xi))

fprintf("%10s | %-9s | %-10s\n","State", 'Iteration', 'Threshold')
for i = 1:n
    fprintf("%10s | %-9s | %-.6f\n", ['DX_', num2str(i)], num2str(numIter(i)), thresholdValue(i));
end

number_label = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

figure('Units', 'centimeters', 'Position', [3 3 17 15]);
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
    set(gca, 'LineWidth',1.2, 'FontSize', fontsize, 'FontName', 'Times');
    ylabel(['X_',num2str(i)], 'FontSize',12);
    
    title(number_label{i}, 'Units', 'normalized', 'Position', [0.04, 0.7, 0], 'FontWeight', 'bold');

    ax = gca;
    grid(ax, 'on');
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.3;
    ax.GridColor = [0.4 0.4 0.4];
end
xlabel('Iterations', 'FontSize',12)
legend(ax1, 'Objective Value', 'Modified Relative Error', 'Location', ...
    'northoutside', 'Orientation', 'horizontal',...
    'Position', [0.25 0.95 0.5 0.05]);
legend(ax1, 'boxoff');

if isOutputToFile
    exportgraphics(t1,[figpath, baseFileName,'_MRE.pdf'],'ContentType','vector', 'BackgroundColor', 'none');
end

coefTrue = cell(1, n);
coefTrue{1} = [-1];
coefTrue{2} = [1];
coefTrue{3} = [-1/3];

xlim_range = [1e-4 5e8];
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

    set(gca, 'LineWidth',1.2, 'FontSize', fontsize, 'FontName', 'Times');

    text(2*xlim_range(1)/10, -i, ['X_',num2str(i)], 'FontSize', fontsize);
    set(gca, 'yTick', [], 'FontName', 'Times');

    if i == 1
        grid minor;
    end
end
xticks([10^(-4), 10^(-2), 10^(0), 10^(2), 10^(4), 10^(6), 10^(8)])
xticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}','10^{6}','10^{8}'})
xlabel('Absolute value of coefficients (logarithmic coordinate)', 'FontSize',fontsize)
legend(gca, 'Initial Least-square Estimation', 'Key Features', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 12);
legend(gca, 'boxoff');
ylim([-(n+1) 0]);
if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_a.pdf'],'ContentType','vector', 'BackgroundColor','none');
end


xlim_range = [1e-5 1e4];
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

    text(3*xlim_range(1)/10, -i, ['X_{',num2str(i),'}'], 'FontSize', fontsize);

    text(lowerBoundary, -i+0.25, str_Pi, 'FontName', 'Times', 'FontWeight', 'Bold', 'FontSize', 12, 'Interpreter','latex');

    set(gca, 'yTick', []);

    if i == 1
        hold on;
        grid minor;
    end
end
xticks([10^(-4), 10^(-3), 10^(-2), 10^(-1), 10^(0), 10^(1), 10^(2), 10^(3)])
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'})
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

text(Pi_min, -n-0.25, str_Pi, 'FontName', 'Times', 'FontWeight', 'Bold', 'FontSize', 12, 'Interpreter','latex');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_b.pdf'],'ContentType','vector', 'BackgroundColor','none');
end

%% integrate true and identified systems
func_name = 'test_3D_Euler_Rigid_Body_Dynamics';
mg = ModelGenerate(Charset, Xi, func_name);

[t_span,x_true]=ode45(@(t,x) erbd(t, x),tspan,x0,options);   % true model
[~, x_identified] = ode45(func_name,tspan,x0,options);  % approximate
delete([func_name,'.m']);

%% Figures
number_label = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

figure('Units', 'centimeters', 'Position', [3 3 17 15]); 
t3 = tiledlayout(n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

set(groot, 'defaultAxesFontName', 'Times New Roman',...  
    'defaultTextFontName', 'Times New Roman',...
    'defaultAxesLabelFontSizeMultiplier', 1.1,...n
    'defaultAxesTitleFontSizeMultiplier', 1.2);

lineStyle = struct(...
    'True', {'Color', [0.2 0.4 0.8], 'LineWidth', 2.5, 'LineStyle', '-'},... 
    'Identified', {'Color', [0.9 0.4 0.1], 'LineWidth', 2, 'LineStyle', '--'}); 

commonYLim = [min(min(x_true)), max(max(x_true))]; 

for i = 1:n
    ax = nexttile;
    plot(t_span, x_true(:,i), lineStyle.True);
    hold on;
    plot(t_span, x_identified(:,i), lineStyle.Identified);
    
    tmp_1 = x_true(:,i);
    tmp_2 = x_identified(:,i);
    fprintf('X_%d: RMSE=%.6f\n', i, error_func(tmp_1, tmp_2, 'rmse'));
    
    tmp_YLim = [min(min(x_true(:,i))), max(max(x_true(:,i)))];
    ylim(yrange_extend(tmp_YLim,0.1,0.05))
    set(ax, 'LineWidth', 1.5,...          
        'FontSize', fontsize,...                
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
        'FontSize', fontsize,...
        'FontWeight', 'normal');

    title(number_label{i},...
        'FontSize', fontsize,...
        'HorizontalAlignment', 'left',...
        'Units', 'normalized',...
        'Position', [0.012 0.88],...       
        'FontWeight', 'bold');
end

xlabel(t3, 'Time (s)',...                 
    'FontSize', fontsize,...
    'FontWeight', 'normal');

leg = legend(ax, {'True','Identified'},...
    'Box', 'off',...
    'Orientation', 'horizontal',...
    'Position', [0.35 0.97 0.5 0.03],...  
    'FontSize', fontsize);
leg.ItemTokenSize = [15,5];               

if isOutputToFile
    exportgraphics(t3,[figpath, baseFileName,'_c.pdf'],'ContentType','vector', 'BackgroundColor', 'none');
end