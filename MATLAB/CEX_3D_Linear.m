clear, clc, close all;
addpath('./utils');
figpath = './figures/';
isOutputToFile = 1; % Switch for outputting figures in PDF format.
baseFileName = 'CEX_3D_Linear';

%% generate Data
polyorder = 2;

n = 3;  % 3D system
A = [-0.0001 0.002 0; -0.002 -0.0001 0 ; 0 0 -0.0003];
rhs = @(x) A*x;   % ODE right hand side
tspan=0:.05:5000;   % time span
x0 = [2; 0; 1];        % initial conditions
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate

% figure;
% tiledlayout(n, 1);
% for i = 1:n
%     nexttile;
%     plot(t, x(:, i))
% end

%% compute Derivative
rng(123456789);
eps = 1e-4;
for i=1:length(x)
    dx_free(i,:) = A*x(i,:)';
end
dx = dx_free + eps*randn(size(dx_free));

figure('Units', 'centimeters', 'Position', [5 5 17 15]);
t0 = tiledlayout(n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:n
    eval(['ax',num2str(i),'=nexttile;']);
    plot(t, dx(:, i))
    hold on;
    plot(t, dx_free(:,i),'--', 'LineWidth', 2)

    set(gca, 'LineWidth',1.5, 'FontSize', 12, 'FontName', 'Times');
    ylabel(['DX_',num2str(i)], 'FontSize',12);
end
xlabel('Time', 'FontSize',12)
legend(ax1, 'Noisy observation', 'Clean observation', 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'FontSize', 12,'Position', [0.25 0.97 0.5 0.02]);
legend(ax1, 'boxoff');

if isOutputToFile
    exportgraphics(t0,[figpath, baseFileName,'_DX.pdf'],'ContentType','vector');
end


%% build library of nonlinear time series
% [Theta, Charset, pad] = getLib(x, polyorder);
% disp(Charset)
% m = size(Theta,2);
[Theta, Charset] = getLibrary(x, polyorder);

Xi_0 = Theta \ dx;
format long;
% tbl_1 = table(Charset, Xi_0);
% writetable(tbl_1);
data_file = 'Charset_Xi.csv';
writetable(table(Charset, Xi_0), data_file);
% delete(data_file); % uncomment this line and run this single line to quick delete the data file.

%% compute Sparse regression: STRidge_ABSS
% [Xi, thresholdValue, numIter, history] = sparsifyDynamics_ABSS(Theta, dx);
[Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);

% disp(thresholdValue)
% disp(numIter)

disp(table(Charset, Xi))

fprintf("%10s | %-9s | %-10s\n","State", 'Iteration', 'Threshold')
for i = 1:n
    fprintf("%10s | %-9s | %-.6f\n", ['DX_', num2str(i)], num2str(numIter(i)), thresholdValue(i));
end

%% comparisions: matrixplot
for k = 1:3
    iter_1 = numIter(k)+1;
    data = history.indmat{k}(:,1:iter_1);
    
    XVarNames = cell(iter_1, 1);
    for i=0:iter_1-1
        XVarNames{i+1} = num2str(i);
    end
    % disp(XVarNames)
    YVarNames = Charset2Tex(Charset);
    matrixplot(data, ...
            'XVarNames',XVarNames, 'YVarNames', YVarNames, 'DisplayOpt','off');
    
    fontname('Times');
    fontsize(15,"points")
    set(gcf, 'Position', [0.72 0.28 0.13 0.2])
    
    if isOutputToFile
        exportgraphics(gcf,[figpath, baseFileName,'_ChangesWithIteration_X',num2str(k),'.pdf'],'ContentType','image');
    end
end

number_label = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

figure('Units', 'centimeters', 'Position', [5 5 17 22]); 
t1 = tiledlayout(n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:n
    eval(['ax',num2str(i),'=nexttile;']);
    yyaxis left;
    y1 = history.objval{i}(:);
    plot(y1, '-', 'LineWidth',2, 'Marker', '*', 'MarkerSize', 10);
    ylim( [yrange_extend(y1)] );
    
    yyaxis right;
    y2 = history.MRE{i}(:);
    plot(y2, '.-', 'LineWidth',2, 'Marker', 'o', 'MarkerSize', 10);
    ylim( [yrange_extend(y2)] );

    grid on;
    set(gca, 'LineWidth',2, 'FontSize', 20, 'FontName', 'Times');
    ylabel(['X_',num2str(i)], 'FontSize',20);

end
xlabel('Iterations', 'FontSize',20)
legend(ax1, 'Objective Value', 'Modified Relative Error', 'Location', 'northoutside', 'Orientation', 'horizontal');
legend(ax1, 'boxoff');

% set(gcf, 'Position', [200, 100, 1200, 900])

if isOutputToFile
    exportgraphics(t1,[figpath, baseFileName,'_MRE.pdf'],'ContentType','vector');
end

coefTrue = cell(1, n);
for i = 1:n
    col = A(i, :);
    nz = col~=0;
    coefTrue{i} = col(nz); 
end

xlim_range = [5e-7 1e-2];
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

    text(4*xlim_range(1)/10, -i, ['X_',num2str(i)]);
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
    
    
    set(gca, 'LineWidth',1.2, 'FontSize', 12, 'FontName', 'Times');

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

dx_identified = Theta * Xi;
figure;
tiledlayout(n, 1);
for i = 1:n
    eval(['ax',num2str(i),'=nexttile;']);
    plot(t, dx(:, i), 'LineWidth',1.5); hold on;
    plot(t, dx_identified(:, i), '--','LineWidth',1.5)

    set(gca, 'LineWidth',1.2, 'FontSize', 12, 'FontName', 'Times');
    ylabel(['DX_',num2str(i)], 'FontSize',12);
end
xlabel('Time', 'FontSize',12)
legend(ax1, 'True', 'Identified', 'Location', 'northoutside', 'Orientation', 'horizontal');
legend(ax1, 'boxoff');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_DX_True_and_Identified.pdf'],'ContentType','vector');
end

%% integrate true and identified systems
% x_true = ode_integration(dx, tspan, x0);
% t_span = tspan;
% [~,x_identified]=ode45(@(t,x)sparseGalerkin_ABSS(t,x,polyorder, Xi),tspan,x0,options);  % approximate
func_name = 'test_CEX_3D_Linear';
mg = ModelGenerate(Charset, Xi, func_name);

[t_span,x_true]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
[~, x_identified] = ode45(func_name,tspan,x0,options);  % approximate
delete([func_name,'.m']);

%% Figures
% number_label = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
% 
% figure;
% t3 = tiledlayout(n, 1);
% t3.TileSpacing = 'compact';
% for i = 1:n
%     eval(['ax',num2str(i),'=nexttile;']);
%     tmp_1 = x_true(:, i);
%     plot(t_span, tmp_1, 'r-', 'LineWidth',1.5);
%     hold on;
%     ylim( [yrange_extend(tmp_1)] );
% 
%     tmp_2 = x_identified(:,i);
%     plot(t_span, tmp_2, 'b--', 'LineWidth',1.5);
%     ylim( [yrange_extend(tmp_2)] );
% 
%     set(gca, 'LineWidth',1.2, 'FontSize', 12, 'FontName', 'Times');
%     ylabel(['X_',num2str(i)], 'FontSize',12);
% 
%     title(number_label{i}, 'Units', 'normalized', 'Position', [0.04, 0.7, 0], 'FontWeight', 'bold');
% 
%     fprintf('X_%d: RMSE=%.6f\n', i, error_func(tmp_1, tmp_2, 'rmse'));
% end
% xlabel('Time', 'FontSize',12)
% legend(ax1, 'True', 'Identified', 'Location', 'northoutside', 'Orientation', 'horizontal');
% legend(ax1, 'boxoff');

%%
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

commonYLim = [min(min(x_true)), max(max(x_true))]; 

for i = 1:n
    ax = nexttile;
    plot(t_span, x_true(:,i), lineStyle.True);
    hold on;
    plot(t_span, x_identified(:,i), lineStyle.Identified);
    
    tmp_1 = x_true(:,i);
    tmp_2 = x_identified(:,i);
    % fprintf('X_%d: RMSE=%.6f\n', i, error_func(tmp_1, tmp_2, 'rmse'));
    % fprintf('X_%d: T=[0, %.2f] RMSE=%.6f\n', i, validPredTime, ...
    %     error_func(tmp_1(1:validPredInd), tmp_2(1:validPredInd), 'rmse'));

    tmp_YLim = [min(min(x_true(:,i))), max(max(x_true(:,i)))];
    ylim(yrange_extend(tmp_YLim,0.1,0.05))
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
        'FontWeight', 'normal')%,...

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

if isOutputToFile
    exportgraphics(t3,[figpath, baseFileName,'_c.pdf'],'ContentType','vector');
end