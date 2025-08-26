clear all, close all, clc
figpath = './figures/';
addpath('./utils');
isOutputToFile = 1; % Switch for outputting figures in PDF format.

%% generate Data
polyorder = 2;

cell_smallest_singular_list = cell(3, 1);
cell_l2_error = cell(3, 1);
cell_fail_SINDy = cell(3, 1);
cell_fail_MRE = cell(3, 1);
cell_fail_FTest = cell(3, 1);
cnt = 0;

MRE_threshold = 0.01;
eta = 1e-8;
title_for_values = 'MRE_001_Reg_minus_8';

dim = 4:23;
for eps = [1, 0.1, 0.01]
    smallest_singular_list = zeros(length(dim), 1);
    l2_error = zeros(length(dim), 1);
    fail_SINDy = zeros(length(dim), 1);
    fail_MRE = zeros(length(dim), 1);
    fail_FTest = zeros(length(dim), 1);
    for k = dim
        n = k;
        F = 8; 
        
        x0 = 8*ones(1, n);
        x0(1) = 1;
        
        % Integrate
        tspan=0:0.01:15;
        N = length(tspan);
        options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
        
        func_name_true = generate_lorenz96(n, F);
        [t,x]=ode45(func_name_true,tspan,x0,options);
        
        % compute Derivative
        rng(123456789);
        dx = zeros(N, n);
        for i=1:length(x)
            dx(i,:) = eval([func_name_true,'(0,x(i,:))']);
        end
        dx = dx + eps*randn(size(dx));
        
        % build library of nonlinear time series
        [Theta, Charset] = getLibrary(x, polyorder);

        Xi_ls = Theta \ dx;
        true_coef = [8 1 1 1]';
    
        for i = 1:k
            sorted_Xi_i = sort(abs(Xi_ls(:,i)), 'des');
            l2_error(k-4+1) = l2_error(k-4+1)+ sqrt(mean(sum((sorted_Xi_i(1:4)-true_coef).^2)));
        end
        l2_error(k-4+1) = l2_error(k-4+1) / k;
        
        % Singular value decomposition
        [U, S, V] = svd(Theta);
        singular_value_list = diag(S);
        smallest_singular_list(k-4+1) = singular_value_list(1)/singular_value_list(end);
        
        % compute Sparse regression: sequential least squares
        lambda = 0.025;      % lambda is our sparsification knob.
        Xi_SINDy = sparsifyDynamics(Theta,dx,lambda,n);
        
        % compute Sparse regression: using MRE as the metric
        % MRE_threshold = 0.1;
        % Xi_MRE = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-8, MRE_threshold);
        
        Xi_MRE = sparsifyDynamics_STRidge_ABSS(Theta, dx, eta, MRE_threshold);
    
        % compute Sparse regression: using FTest as the metric
        Xi_FTest = sparsifyDynamics_ABSS_FTest(Theta, dx);
    
        fail_SINDy(k-4+1) = sum(sum(Xi_SINDy~=0)~=4)/n;
        fail_MRE(k-4+1) = sum(sum(Xi_MRE~=0)~=4)/n;
        fail_FTest(k-4+1) = sum(sum(Xi_FTest~=0)~=4)/n;
    end

        cnt = cnt + 1;
        cell_smallest_singular_list{cnt} = smallest_singular_list;
        cell_l2_error{cnt} = l2_error;
        cell_fail_SINDy{cnt} = fail_SINDy;
        cell_fail_MRE{cnt} = fail_MRE;
        cell_fail_FTest{cnt} = fail_FTest;
        % 
end

delete('lorenz96_n*');
%% 
font_size = 12;
number_label = {'Ⅰ', 'Ⅱ', 'Ⅲ'};
figure;
t1 = tiledlayout(3, 1);
t1.TileSpacing = 'compact';

ax1 = nexttile;
ind = 1;
plot(dim, 1-cell_fail_SINDy{ind},'o-', 'LineWidth', 2, 'MarkerSize', 5);
hold on;
plot(dim, 1-cell_fail_MRE{ind},'s--', 'LineWidth', 2, 'MarkerSize', 5);
plot(dim, 1-cell_fail_FTest{ind},'^:', 'LineWidth', 2, 'MarkerSize', 5);

ylabel('$S$', 'FontSize', 12, 'Interpreter', 'latex');

set(gca, 'LineWidth', 1.5, 'FontSize', font_size, 'FontName', 'Times');
title(number_label{ind}, 'Units', 'normalized', 'Position', [0.04, 0.65, 0], 'FontWeight', 'bold');

tmp = cell_fail_SINDy{ind};
sr = find(tmp==0);
x_pos = 3 + sr(end);
line([x_pos, x_pos], [-0.1, 1.1], 'lineStyle', '--', 'LineWidth', 2, 'HandleVisibility','off');

xlim([3 dim(end)+1]);
xticks(3:1:dim(end)+1);
ylim([-0.1 1.1]);

ax2 = nexttile;
ind = 2;
plot(dim, 1-cell_fail_SINDy{ind},'o-', 'LineWidth', 2, 'MarkerSize', 5);
hold on;
plot(dim, 1-cell_fail_MRE{ind},'s--', 'LineWidth', 2, 'MarkerSize', 5);
plot(dim, 1-cell_fail_FTest{ind},'^:', 'LineWidth', 2, 'MarkerSize', 5);

ylabel('$S$', 'FontSize', 12, 'Interpreter', 'latex');

set(gca, 'LineWidth', 1.5, 'FontSize', font_size, 'FontName', 'Times');
title(number_label{ind}, 'Units', 'normalized', 'Position', [0.04, 0.65, 0], 'FontWeight', 'bold');

tmp = cell_fail_SINDy{ind};
sr = find(tmp==0);
x_pos = 3 + sr(end);
line([x_pos, x_pos], [-0.1, 1.1], 'lineStyle', '--', 'LineWidth', 2, 'HandleVisibility','off');

xlim([3 dim(end)+1]);
xticks(3:1:dim(end)+1);
ylim([-0.1 1.1]);

ax3 = nexttile;
ind = 3;
plot(dim, 1-cell_fail_SINDy{ind},'o-', 'LineWidth', 2, 'MarkerSize', 5);
hold on;
plot(dim, 1-cell_fail_MRE{ind},'s--', 'LineWidth', 2, 'MarkerSize', 5);
plot(dim, 1-cell_fail_FTest{ind},'^:', 'LineWidth', 2, 'MarkerSize', 5);

xlabel('$n$', 'FontSize', 12, 'Interpreter', 'latex')
ylabel('$S$', 'FontSize', 12, 'Interpreter', 'latex');

set(gca, 'LineWidth', 1.5, 'FontSize', font_size, 'FontName', 'Times');
title(number_label{ind}, 'Units', 'normalized', 'Position', [0.04, 0.65, 0], 'FontWeight', 'bold');

tmp = cell_fail_SINDy{ind};
sr = find(tmp==0);
x_pos = 3 + sr(end);
line([x_pos, x_pos], [-0.1, 1.1], 'lineStyle', '--', 'LineWidth', 2, 'HandleVisibility','off');

xlim([3 dim(end)+1]);
xticks(3:1:dim(end)+1);
ylim([-0.1 1.1]);

legend(ax1, 'SINDy', 'BSS with MRE', 'BSS with FTest', 'Location', 'northoutside', 'Orientation', 'horizontal');
legend(ax1, 'boxoff');


if isOutputToFile
    exportgraphics(gcf,[figpath, 'EX_HD_Lorenz_96_', title_for_values, '_SRR.pdf'],'ContentType','vector', 'BackgroundColor','none');
end
