clear, clc, close all;
figpath = './figures/';
isOutputToFile = 0; % Switch for outputting figures in PDF format.
baseFileName = 'TimeConsumptionAnalysis';

%%
VariableNames = {'Model', 'Mean', 'Std', 'MeanPerStep', 'StdPerStep'};

load TimeConsumptionAnalysisSRBased_Low.mat
T_SRBasedLD = table(model, record(:, 1), record(:, 2), record(:, 3), record(:, 4), 'VariableNames', VariableNames);

load TimeConsumptionAnalysisFSINDyBased_Low.mat
T_FSINDyEnhancedLD = table(model, record(:, 1), record(:, 2), record(:, 3), record(:, 4), 'VariableNames', VariableNames);

load TimeConsumptionAnalysisSRBased_High_Lorenz96.mat
T_SRBasedHighLorenz96 = table(model, record(:, 1), record(:, 2), record(:, 3), record(:, 4), 'VariableNames', VariableNames);

load TimeConsumptionAnalysisFSINDyBased_High_Lorenz96.mat
T_FSINDyEnhancedHighLorenz96 = table(model, record(:, 1), record(:, 2), record(:, 3), record(:, 4), 'VariableNames', VariableNames);

load TimeConsumptionAnalysisSRBased_High.mat
T_SRBasedHD = table(model, record(:, 1), record(:, 2), record(:, 3), record(:, 4), 'VariableNames', VariableNames);

load TimeConsumptionAnalysisFSINDyBased_High.mat
T_FSINDyEnhancedHD = table(model, record(:, 1), record(:, 2), record(:, 3), record(:, 4), 'VariableNames', VariableNames);


dt = 0.01;
T_list = [25, 25, 25, 50, 50, 50]'/dt;
Dim = [2, 2, 2, 3, 3, 3]';

% %% =================================================================================================
figure('Units', 'centimeters', 'Position', [2 2 17 17]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

FontSize = 12;
InnerLineWidth = 1.5;
OuterLineWidth = 1.2;

%% subfig: (a)
nexttile;
bar(1:length(T_SRBasedLD.Model), T_SRBasedLD.Mean, 'EdgeColor', 'none', 'FaceColor', '#6888F5');
hold on;
er = errorbar(T_SRBasedLD.Mean, T_SRBasedLD.Std, 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(T_SRBasedLD.Model)
ylabel('Average simulation time (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(a)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

%% subfig: (b)
nexttile;
bar(1:length(T_SRBasedLD.Model), T_SRBasedLD.Mean./(T_list.*Dim), 'EdgeColor', 'none', 'FaceColor', '#6888F5');

hold on;
er = errorbar(T_SRBasedLD.Mean./(T_list.*Dim), T_SRBasedLD.Std./(T_list.*Dim), 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(T_SRBasedLD.Model)
ylabel('Average simulation time per step (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(b)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

%% subfig: (c)
nexttile;
b2 = bar(1:length(T_FSINDyEnhancedLD.Model), T_FSINDyEnhancedLD.Mean, 'EdgeColor', 'none', 'FaceColor', '#D77071');
hold on;
er = errorbar(T_FSINDyEnhancedLD.Mean, T_FSINDyEnhancedLD.Std, 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(T_FSINDyEnhancedLD.Model)
ylabel('Average simulation time (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(c)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');
%% subfig: (d)
nexttile;
bar(1:length(T_FSINDyEnhancedLD.Model), T_FSINDyEnhancedLD.Mean./(T_list.*Dim), 'EdgeColor', 'none', 'FaceColor', '#D77071');

hold on;
er = errorbar(T_FSINDyEnhancedLD.Mean./(T_list.*Dim), T_FSINDyEnhancedLD.Std./(T_list.*Dim), 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(T_FSINDyEnhancedLD.Model)
ylabel('Average simulation time per step (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(d)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_Summary_LowDimension.pdf'],'ContentType','vector', 'BackgroundColor', 'none');
end

%%
dt = 0.01;
T_list = [500, 200, 200, 15]'/dt;
Dim = [4, 4, 6, 8]';

figure('Units', 'centimeters', 'Position', [2 2 17 17]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

FontSize = 11;
InnerLineWidth = 1.5;
OuterLineWidth = 1.2;

%% subfig: (a)
nexttile;
bar(1:length(T_SRBasedHD.Model), T_SRBasedHD.Mean, 'EdgeColor', 'none', 'FaceColor', '#6888F5');
hold on;
er = errorbar(T_SRBasedHD.Mean, T_SRBasedHD.Std, 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(T_SRBasedHD.Model)
ylabel('Average simulation time (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(a)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

%% subfig: (b)
nexttile;
bar(1:length(T_SRBasedHD.Model), T_SRBasedHD.Mean./(T_list.*Dim), 'EdgeColor', 'none', 'FaceColor', '#6888F5');

hold on;
er = errorbar(T_SRBasedHD.Mean./(T_list.*Dim), T_SRBasedHD.Std./(T_list.*Dim), 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(T_SRBasedHD.Model)
ylabel('Average simulation time per step (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(b)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

%% subfig: (c)
nexttile;
bar(1:length(T_FSINDyEnhancedHD.Model), T_FSINDyEnhancedHD.Mean, 'EdgeColor', 'none', 'FaceColor', '#D77071');

hold on;
er = errorbar(1:length(T_FSINDyEnhancedHD.Model), T_FSINDyEnhancedHD.Mean, T_FSINDyEnhancedHD.Std, 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(T_FSINDyEnhancedHD.Model)
ylabel('Average simulation time (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(c)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

%% subfig: (d)
nexttile;
bar(1:length(T_FSINDyEnhancedHD.Model), T_FSINDyEnhancedHD.Mean./(T_list.*Dim), 'EdgeColor', 'none', 'FaceColor', '#D77071');

hold on;
er = errorbar(1:length(T_FSINDyEnhancedHD.Model), T_FSINDyEnhancedHD.Mean./(T_list.*Dim), T_FSINDyEnhancedHD.Std./(T_list.*Dim), 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(T_FSINDyEnhancedHD.Model)
ylabel('Average simulation time per step (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(d)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_Summary_HighDimension.pdf'],'ContentType','vector', 'BackgroundColor', 'none');
end

%%
dt = 0.01;
T_list = 15*ones(17,1)/dt;
Dim = [4:20]';

figure('Units', 'centimeters', 'Position', [2 2 17 17]);
tiledlayout(2,4,'TileSpacing', 'compact','Padding','compact');

FontSize = 12;
InnerLineWidth = 1.5;
OuterLineWidth = 1.2;

%% subfig: (a)
nexttile([1 2]);
bar(4:20, T_SRBasedHighLorenz96.Mean, 'EdgeColor', 'none', 'FaceColor', '#6888F5');

hold on;
er = errorbar(4:20, T_SRBasedHighLorenz96.Mean, T_SRBasedHighLorenz96.Std, 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticks(4:20)
xlabel('Dimension')
ylabel('Average simulation time (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(a)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

%% subfig: (b)
nexttile([1 2]);
bar(4:20, T_SRBasedHighLorenz96.Mean./(T_list.*Dim), 'EdgeColor', 'none', 'FaceColor', '#6888F5');

hold on;
er = errorbar(4:20, T_SRBasedHighLorenz96.Mean./(T_list.*Dim), T_SRBasedHighLorenz96.Std./(T_list.*Dim), 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticks(4:20)
xlabel('Dimension')
ylabel('Average simulation time per step (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(b)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

%% subfig: (c)
nexttile([1 2]);
bar(4:20, T_FSINDyEnhancedHighLorenz96.Mean, 'EdgeColor', 'none', 'FaceColor', '#D77071');

hold on;
er = errorbar(4:20, T_FSINDyEnhancedHighLorenz96.Mean, T_FSINDyEnhancedHighLorenz96.Std, 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticks(4:20)
xlabel('Dimension')
ylabel('Average simulation time (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)

title('(c)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

%% subfig: (d)
nexttile([1 2]);
bar(4:20, T_FSINDyEnhancedHighLorenz96.Mean./(T_list.*Dim), 'EdgeColor', 'none', 'FaceColor', '#D77071');

hold on;
er = errorbar(4:20, T_FSINDyEnhancedHighLorenz96.Mean./(T_list.*Dim), T_FSINDyEnhancedHighLorenz96.Std./(T_list.*Dim), 'LineWidth', InnerLineWidth);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticks(4:20)
ylim([0 1.2e-5]);
xlabel('Dimension')
ylabel('Average simulation time per step (s)')
set(gca, 'FontSize', FontSize, 'FontName', 'Times', 'LineWidth', OuterLineWidth)


title('(d)', 'Units', 'normalized', 'Position', [0.055, 0.91, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_Summary_HighLorenz96.pdf'],'ContentType','vector', 'BackgroundColor', 'none');
end

%% 
figure('Units', 'centimeters', 'Position', [2 2 17 17]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

FontSize = 12;
InnerLineWidth = 1.5;
OuterLineWidth = 1.2; 

nexttile;
semilogy(1:length(T_SRBasedLD.Model), T_SRBasedLD.Mean ./ T_FSINDyEnhancedLD.Mean, ...
    '-o', 'LineWidth', InnerLineWidth, 'MarkerSize', 10, 'MarkerFaceColor','#D77071','MarkerEdgeColor','none');


box off
grid on

xticks(1:6);
xlim([0 7]);
ylim([1e1 1e3]);
xticklabels(T_SRBasedLD.Model);

ylabel('Speed-up ratio')

set(gca, 'LineWidth', 1.5,...          
    'FontSize', FontSize,...                
    'FontName', 'Times',...
    'TickDir', 'out',...              
    'TickLength', [0.015 0.015],...   
    'YMinorTick', 'on',...
    'Box', 'off');                    

title('(a)', 'Units', 'normalized', 'Position', [0.075, 0.93, 0], 'FontWeight', 'bold');

InnerLineWidth = 1.5;
OuterLineWidth = 1.2; 

nexttile;
semilogy(1:length(T_SRBasedHD.Model), T_SRBasedHD.Mean ./ T_FSINDyEnhancedHD.Mean,...
    '-o', 'LineWidth', InnerLineWidth, 'MarkerSize', 10, 'MarkerFaceColor','#D77071','MarkerEdgeColor','none');

box off
grid on

xticks(1:4);
xlim([0 5]);
ylim([1e1 1e3]);
xticklabels(T_SRBasedHD.Model);

ylabel('Speed-up ratio')

set(gca, 'LineWidth', 1.5,...          
    'FontSize', FontSize,...                
    'FontName', 'Times',...
    'TickDir', 'out',...              
    'TickLength', [0.015 0.015],...   
    'YMinorTick', 'on',...
    'Box', 'off');                    

title('(b)', 'Units', 'normalized', 'Position', [0.075, 0.93, 0], 'FontWeight', 'bold');

InnerLineWidth = 1.5;
OuterLineWidth = 1.2; 

nexttile([1 2]);
semilogy(4:20, T_SRBasedHighLorenz96.Mean ./ T_FSINDyEnhancedHighLorenz96.Mean, ...
    '-o', 'LineWidth', InnerLineWidth, 'MarkerSize', 10, 'MarkerFaceColor','#D77071','MarkerEdgeColor','none');

box off
grid on
xlim([3 21])
xticks(4:20)
ylim([1e2 5e3]);

xlabel('Dimension')
ylabel('Speed-up ratio')

set(gca, 'LineWidth', 1.5,...          
    'FontSize', FontSize,...                
    'FontName', 'Times',...
    'TickDir', 'out',...              
    'TickLength', [0.015 0.015],...   
    'YMinorTick', 'on',...
    'Box', 'off');                    

title('(c)', 'Units', 'normalized', 'Position', [0.035, 0.93, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_Summary_SpeedUpRatio.pdf'],'ContentType','vector', 'BackgroundColor', 'none');
end