clear, clc, close all;
addpath('./utils');

figpath = './figures/';
isOutputToFile = 1;
baseFileName = 'TimeConsumptionAnalysisSRBased';
isOverwrite = 1;

if ~exist("TimeConsumptionAnalysisSRBased_High_Lorenz96.mat") || isOverwrite
    num_of_simulation = 20;
    model = cell(20-4+1, 1);
    record = zeros(20-4+1, 4);
    
    for n = 4:20
        clear x dx
        model_name = ['Lorenz96_',num2str(n),'D'];
        time_record = zeros(num_of_simulation, 1);
        time_record_per_step = zeros(num_of_simulation, 1);
        for k = 1:num_of_simulation
            % generate Data
            polyorder = 2;
            usesine = 0;
            F = 8;
            
            x0 = 8*ones(1, n);
            x0(1) = 1;
            
            % x0=[-8; 8; 27];  % Initial condition
            
            % Integrate
            tspan=[0:0.01:15];
            N = length(tspan);
            options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
            % [t,x]=ode45(@(t,x) lorenz_96(t,x),tspan,x0,options);
            
            
            % [t,x]=ode45(@lorenz96_m_4_F_8,tspan,x0,options);
            
            func_name_true = generate_lorenz96(n, F);
            
            [t,x]=ode45(func_name_true,tspan,x0,options);
            
            % compute Derivative
            rng(123456789);
            dx = zeros(N, n);
            eps = 0.001;
            % eps = 0;
            for i=1:length(x)
                dx(i,:) = eval([func_name_true,'(0,x(i,:))']);
            end
            dx = dx + eps*randn(size(dx));
            
            % build library of nonlinear time series
            [Theta, Charset] = getLibrary(x, polyorder);
            % disp(Charset)
            m = size(Theta,2);
            
            [Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);
            % integrate true and identified systems
            % func_name = ['func_TimeConsumptionAnalysis_',model_name];
            % mg = ModelGenerate(Charset, Xi, func_name);
            % disp(mg);
            
            %[t_span,x_true]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
            tic; 
            time_start = toc; 
            % [~,x_identified]=ode45(func_name,tspan,x0,options);  % approximate
            [~,x_identified]=ode45(@(t,x)sparseGalerkin_ABSS(t,x,polyorder, Xi),tspan,x0,options);  % approximate
            time_end = toc;
            time_elapsed = time_end - time_start;
            time_record(k) = time_elapsed;
            time_record_per_step(k) = time_elapsed/length(tspan);
        end
        mean_time = mean(time_record);
        std_time = std(time_record);
        mean_time_per_step = mean(time_record_per_step);
        std_time_per_step = std(time_record_per_step);
        model{n-4+1} = [num2str(n),'D'];
        record(n-4+1, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];

        fprintf('n = %d, finished!\n', n);
    end
    save('TimeConsumptionAnalysisSRBased_High_Lorenz96.mat','model', 'record');
else
    load TimeConsumptionAnalysisSRBased_High_Lorenz96.mat
end

%% Final Analysis
VariableNames = {'Model', 'Mean', 'Std', 'MeanPerStep', 'StdPerStep'};
T = table(model, record(:, 1), record(:, 2), record(:, 3), record(:, 4), 'VariableNames', VariableNames);

figure;
bar(4:20, T.Mean, 'EdgeColor', 'none', 'FaceColor', '#6888F5');

hold on;
er = errorbar(4:20, T.Mean, T.Std, 'LineWidth', 1.5);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on

xticks(4:20)
xlabel('Dimension')
ylabel('Average simulation time (s)')
set(gca, 'FontSize', 12, 'FontName', 'Times', 'LineWidth', 1)

title('(c)', 'Units', 'normalized', 'Position', [0.04, 0.93, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_c.pdf'],'ContentType','vector');
end

figure;
bar(4:20, T.MeanPerStep, 'EdgeColor', 'none', 'FaceColor', '#6888F5');

hold on;
er = errorbar(4:20, T.MeanPerStep, T.StdPerStep, 'LineWidth', 1.5);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticks(4:20)
xlabel('Dimension')
ylabel('Average simulation time per step (s)')
set(gca, 'FontSize', 12, 'FontName', 'Times', 'LineWidth', 1)

title('(d)', 'Units', 'normalized', 'Position', [0.04, 0.93, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_d.pdf'],'ContentType','vector');
end
