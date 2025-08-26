clear, clc, close all;
addpath('./utils');

figpath = './figures/';
isOutputToFile = 1;
baseFileName = 'TimeConsumptionAnalysisSRBased_High';
isOverwrite = 1;

if ~exist("TimeConsumptionAnalysisSRBased_High.mat") || isOverwrite
    num_of_simulation = 20;
    model = cell(4, 1);
    record = zeros(4, 4);
    
    cnt = 0;
    %% 1
    clear x dx_free dx
    cnt = cnt + 1;
    model_name = '4D_LV';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 2;  % search space up to fifth order polynomials
        usesine = 0;    % no trig functions
        n = 4;

        x0=[0.5; 0.2; 0.1; 0.5];  % Initial condition

        % Integrate
        tspan=[0:0.01:500];
        N = length(tspan);
        options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
        [t,x]=ode45(@(t,x) Hyper_LV(t,x),tspan,x0,options);
        
        % compute Derivative
        rng(123456789);
        eps = 0.001;
        for i=1:length(x)
            dx_free(i,:) = Hyper_LV(0,x(i,:));
        end
        % dx_noise = dx + eps*randn(size(dx));
        dx = dx_free + eps*randn(size(dx_free));
        
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
    model{cnt} = 'Lotka–Volterra (4D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    
    %% 2
    clear x dx_free dx
    cnt = cnt + 1;
    model_name = '4D_Rossler';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 2;  % search space up to fifth order polynomials
        usesine = 0;    % no trig functions
        a = 0.25;
        b = 3;
        c = 0.05;
        d = 0.5;
        
        n = 4;
        
        x0=[-6; 0; 0.5; 14];  % Initial condition
        
        % Integrate
        tspan=[0:0.01:200];
        N = length(tspan);
        options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
        [t,x]=ode45(@(t,x) Hyper_Rossler(t,x,a,b,c,d),tspan,x0,options);
        
        %% compute Derivative
        rng(123456789);
        eps = 0.005;
        for i=1:length(x)
            dx_free(i,:) = Hyper_Rossler(0,x(i,:),a,b,c,d);
        end
        % dx_noise = dx + eps*randn(size(dx));
        dx = dx_free + eps*randn(size(dx_free));
        
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
    model{cnt} = 'Rossler Hyperchaotic (4D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    %% 3
    clear x dx_free dx
    cnt = cnt + 1;
    model_name = '6D_Coupled_Lorenz';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 2;  % search space up to fifth order polynomials
        usesine = 0;    % no trig functions
        n = 6;          % 2D system
        K = 0.5;
        R = 1;
        x0=[0; 0.8; 0.4; 0; 1; 0];  % Initial condition

        % Integrate
        tspan=[0:0.01:200];
        N = length(tspan);
        options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
        [t,x]=ode45(@(t,x) Hyper_Coupled_Lorenz(t,x,K,R),tspan,x0,options);
        
        %% compute Derivative
        rng(123456789);
        eps = 0.005;
        for i=1:length(x)
            dx_free(i,:) = Hyper_Coupled_Lorenz(0,x(i,:),K,R);
        end
        % dx_noise = dx + eps*randn(size(dx));
        dx = dx_free + eps*randn(size(dx_free));
        
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
    model{cnt} = 'Coupled Lorenz (6D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    
    %% 4
    clear x dx_free dx
    cnt = cnt + 1;
    model_name = '8D_Lorenz96';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 2;
        usesine = 0;
        F = 8;
        n = 8;
            
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
        eps = 0.01;
        % eps = 0;
        for i=1:length(x)
            dx_free(i,:) = eval([func_name_true,'(0,x(i,:))']);
        end
        dx = dx_free + eps*randn(size(dx_free));
        
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
    model{cnt} = 'Lorenz96 (8D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    save('TimeConsumptionAnalysisSRBased_High.mat','model', 'record');
else
    load TimeConsumptionAnalysisSRBased_High.mat
end

%% Final Analysis
dt = 0.01;
T_list = [500, 200, 200, 15]'/dt;
Nf = [15, 15, 28, 45]';
Dim = [4, 4, 6, 8]';

VariableNames = {'Model', 'Mean', 'Std', 'MeanPerStep', 'StdPerStep'};
T = table(model, record(:, 1), record(:, 2), record(:, 3), record(:, 4), 'VariableNames', VariableNames);

figure;
bar(1:length(model), T.Mean, 'EdgeColor', 'none', 'FaceColor', '#6888F5');

hold on;
er = errorbar(T.Mean, T.Std, 'LineWidth', 1.5);
er.Color = [0 0 0];
er.LineStyle = 'none';
  
box off
grid on
xticklabels(model)
ylabel('Average simulation time (s)')
% set(gca, 'TickLabelInterpreter', 'none');
set(gca, 'FontSize', 12, 'FontName', 'Times', 'LineWidth', 1)
% set(gca, 'YScale', 'log');

title('(a)', 'Units', 'normalized', 'Position', [0.04, 0.93, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_a.pdf'],'ContentType','vector');
end

figure;
bar(1:length(model), T.MeanPerStep./(T_list.*Dim), 'EdgeColor', 'none', 'FaceColor', '#6888F5');

hold on;
er = errorbar(T.MeanPerStep./(T_list.*Dim), T.StdPerStep./(T_list.*Dim), 'LineWidth', 1.5);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(model)
ylabel('Average simulation time per step (s)')
% set(gca, 'TickLabelInterpreter', 'none');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 12, 'FontName', 'Times', 'LineWidth', 1)

title('(b)', 'Units', 'normalized', 'Position', [0.04, 0.93, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_b.pdf'],'ContentType','vector');
end
