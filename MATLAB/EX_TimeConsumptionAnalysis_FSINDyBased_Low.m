clear, clc, close all;
addpath('./utils');

figpath = './figures/';
isOutputToFile = 1;
baseFileName = 'TimeConsumptionAnalysis';
isOverwrite = 1;

if ~exist("TimeConsumptionAnalysisFSINDyBased_Low.mat") || isOverwrite
    num_of_simulation = 200;
    model = cell(6, 1);
    record = zeros(6, 4);
    
    cnt = 0;
    %% 1
    clear x dx
    cnt = cnt + 1;
    model_name = '2D_Linear';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 3;  % search space up to fifth order polynomials
        usesine = 0;    % no trig functions
        n = 2;          % 2D system
        A = [-.1 2; -2 -.1];  % dynamics
        rhs = @(x)A*x;   % ODE right hand side
        tspan=[0:.01:25];   % time span
        x0 = [2; 0];        % initial conditions
        options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
        [t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate
        
        % compute Derivative
        rng(123456789);
        eps = .05;      % noise strength
        for i=1:length(x)
            dx(i,:) = A*x(i,:)';
        end
        dx = dx + eps*randn(size(dx));   % add noise
        
        % build library of nonlinear time series
        % [Theta, Charset, pad] = getLib(x, polyorder);
        [Theta, Charset] = getLibrary(x, polyorder);
        % disp(Charset)
        m = size(Theta,2);
        
        [Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);
        % integrate true and identified systems
        func_name = ['func_TimeConsumptionAnalysis_',model_name];
        mg = ModelGenerate(Charset, Xi, func_name);
        % disp(mg);
        
        %[t_span,x_true]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
        tic; 
        time_start = toc; 
        [~,x_identified]=ode45(func_name,tspan,x0,options);  % approximate
        time_end = toc;
        delete([func_name,'.m']);
        time_elapsed = time_end - time_start;
        time_record(k) = time_elapsed;
        time_record_per_step(k) = time_elapsed/length(tspan);
    end
    mean_time = mean(time_record);
    std_time = std(time_record);
    mean_time_per_step = mean(time_record_per_step);
    std_time_per_step = std(time_record_per_step);
    model{cnt} = 'Linear (2D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    
    %% 2
    clear x dx
    cnt = cnt + 1;
    model_name = '2D_Cubic';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 5;  % search space up to fifth order polynomials
        usesine = 0;    % no trig functions
        n = 2;          % 2D system
        A = [-.1 2; -2 -.1];
        rhs = @(y)A*y.^3;   % ODE right hand side
        tspan=[0:.01:25];   % time span
        x0 = [2; 0];        % initial conditions
        options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
        [t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate
        
        % compute Derivative
        rng(123456789);
        eps = .05;      % noise strength
        for i=1:length(x)
            dx(i,:) = A*(x(i,:).^3)';
        end
        dx = dx + eps*randn(size(dx));   % add noise
        
        % build library of nonlinear time series
        [Theta, Charset, pad] = getLib(x, polyorder);
        % disp(Charset)
        m = size(Theta,2);
        
        [Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);
        % integrate true and identified systems
        func_name = ['func_TimeConsumptionAnalysis_',model_name];
        mg = ModelGenerate(Charset, Xi, func_name);
        % disp(mg);
        
        %[t_span,x_true]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
        tic; 
        time_start = toc; 
        [~,x_identified]=ode45(func_name,tspan,x0,options);  % approximate
        time_end = toc;
        time_elapsed = time_end - time_start;
        time_record(k) = time_elapsed;
        time_record_per_step(k) = time_elapsed/length(tspan);
    end
    mean_time = mean(time_record);
    std_time = std(time_record);
    mean_time_per_step = mean(time_record_per_step);
    std_time_per_step = std(time_record_per_step);
    model{cnt} = 'Cubic (2D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    %% 3
    clear x dx
    cnt = cnt + 1;
    model_name = '2D_vanDerPol';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 3;  % search space up to fifth order polynomials
        usesine = 0;    % no trig functions
        n = 2;          % 2D system
        
        u = 2;
        rhs = @(x) [x(2);
                    -x(1) - u * (x(1)^2 - 1) * x(2)];
        % X0 = [0.1 0];
        x0 = [0.01 0];
        
        % A = [-.1 2; -2 -.1];  % dynamics
        % rhs = @(x)A*x;   % ODE right hand side
        tspan=[0:0.01:25];   % time span
        % x0 = [2; 0];        % initial conditions
        options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
        [t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate
        
        % compute Derivative
        rng(123456789);
        eps = .05;      % noise strength
        for i=1:length(x)
            dx(i,:) = rhs(x(i,:))';
        end
        dx = dx + eps*randn(size(dx));   % add noise
        
        % build library of nonlinear time series
        [Theta, Charset, pad] = getLib(x, polyorder);
        % disp(Charset)
        m = size(Theta,2);
        
        [Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);
        % integrate true and identified systems
        func_name = ['func_TimeConsumptionAnalysis_',model_name];
        mg = ModelGenerate(Charset, Xi, func_name);
        % disp(mg);
        
        %[t_span,x_true]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
        tic; 
        time_start = toc; 
        [~,x_identified]=ode45(func_name,tspan,x0,options);  % approximate
        time_end = toc;
        time_elapsed = time_end - time_start;
        time_record(k) = time_elapsed;
        time_record_per_step(k) = time_elapsed/length(tspan);
    end
    mean_time = mean(time_record);
    std_time = std(time_record);
    mean_time_per_step = mean(time_record_per_step);
    std_time_per_step = std(time_record_per_step);
    model{cnt} = 'van Der Pol (2D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    
    %% 4
    clear x dx
    cnt = cnt + 1;
    model_name = '3D_Linear';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 2;
        usesine = 0;
        n = 3;  % 3D system
        A = [-.1 2 0; -2 -.1 0 ; 0 0 -.3];
        rhs = @(x)A*x;   % ODE right hand side
        tspan=[0:.01:50];   % time span
        x0 = [2; 0; 1];        % initial conditions
        options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
        [t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate
    
        % compute Derivative
        rng(123456789);
        eps = 0.01;
        for i=1:length(x)
            dx(i,:) = A*x(i,:)';
        end
        dx = dx + eps*randn(size(dx));
        
        % build library of nonlinear time series
        [Theta, Charset, pad] = getLib(x, polyorder);
        % disp(Charset)
        m = size(Theta,2);
        
        [Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);
        % integrate true and identified systems
        func_name = ['func_TimeConsumptionAnalysis_',model_name];
        mg = ModelGenerate(Charset, Xi, func_name);
        % disp(mg);
        
        %[t_span,x_true]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
        tic; 
        time_start = toc; 
        [~,x_identified]=ode45(func_name,tspan,x0,options);  % approximate
        time_end = toc;
        time_elapsed = time_end - time_start;
        time_record(k) = time_elapsed;
        time_record_per_step(k) = time_elapsed/length(tspan);
    end
    mean_time = mean(time_record);
    std_time = std(time_record);
    mean_time_per_step = mean(time_record_per_step);
    std_time_per_step = std(time_record_per_step);
    model{cnt} = 'Linear (3D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    %% 5
    clear x dx
    cnt = cnt + 1;
    model_name = '3D_Rossler';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 3;
        usesine = 0;
        
        alpha = 0.2;  % Rossler's parameters (chaotic)
        beta = 5.7;
        
        n = 3;
        
        x0=[2.9; -1.3; 25];  % Initial condition
        
        % Integrate
        tspan=[0:0.01:50];
        N = length(tspan);
        options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
        [t,x]=ode45(@(t,x) rossler(t,x,alpha,beta),tspan,x0,options);
    
    
        % compute Derivative
        rng(123456789);
        eps = 1;
        for i=1:length(x)
            dx(i,:) = rossler(0,x(i,:),alpha,beta);
        end
        dx = dx + eps*randn(size(dx));
        
        % build library of nonlinear time series
        [Theta, Charset, pad] = getLib(x, polyorder);
        % disp(Charset)
        m = size(Theta,2);
        
        [Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);
        % integrate true and identified systems
        func_name = ['func_TimeConsumptionAnalysis_',model_name];
        mg = ModelGenerate(Charset, Xi, func_name);
        % disp(mg);
        
        %[t_span,x_true]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
        tic; 
        time_start = toc; 
        [~,x_identified]=ode45(func_name,tspan,x0,options);  % approximate
        time_end = toc;
        time_elapsed = time_end - time_start;
        time_record(k) = time_elapsed;
        time_record_per_step(k) = time_elapsed/length(tspan);
    end
    mean_time = mean(time_record);
    std_time = std(time_record);
    mean_time_per_step = mean(time_record_per_step);
    std_time_per_step = std(time_record_per_step);
    model{cnt} = 'Rossler (3D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    %% 6
    clear x dx
    cnt = cnt + 1;
    model_name = '3D_Lorenz63';
    time_record = zeros(num_of_simulation, 1);
    time_record_per_step = zeros(num_of_simulation, 1);
    for k = 1:num_of_simulation
        % generate Data
        polyorder = 3;
        usesine = 0;
        
        sigma = 10;  % Lorenz's parameters (chaotic)
        beta = 8/3;
        rho = 28;
        
        n = 3;
        
        x0=[-8; 8; 27];  % Initial condition
        
        % Integrate
        tspan=[0:0.01:50];
        N = length(tspan);
        options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
        [t,x]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options);
    
    
        % compute Derivative
        rng(123456789);
        eps = 1;
        for i=1:length(x)
            dx(i,:) = lorenz(0,x(i,:),sigma,beta,rho);
        end
        dx = dx + eps*randn(size(dx));
        
        % build library of nonlinear time series
        [Theta, Charset, pad] = getLib(x, polyorder);
        % disp(Charset)
        m = size(Theta,2);
        
        [Xi, thresholdValue, numIter, history] = sparsifyDynamics_STRidge_ABSS(Theta, dx, 1e-6);
        % integrate true and identified systems
        func_name = ['func_TimeConsumptionAnalysis_',model_name];
        mg = ModelGenerate(Charset, Xi, func_name);
        % disp(mg);
        
        %[t_span,x_true]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
        tic; 
        time_start = toc; 
        [~,x_identified]=ode45(func_name,tspan,x0,options);  % approximate
        time_end = toc;
        time_elapsed = time_end - time_start;
        time_record(k) = time_elapsed;
        time_record_per_step(k) = time_elapsed/length(tspan);
    end
    mean_time = mean(time_record);
    std_time = std(time_record);
    mean_time_per_step = mean(time_record_per_step);
    std_time_per_step = std(time_record_per_step);
    model{cnt} = 'Lorenz 63 (3D)';
    record(cnt, :) = [mean_time, std_time, mean_time_per_step, std_time_per_step];
    
    save('TimeConsumptionAnalysisFSINDyBased_Low.mat','model', 'record');
else
    load TimeConsumptionAnalysisFSINDyBased_Low.mat
end

delete('func_TimeConsumptionAnalysis_*');

%% Final Analysis
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
set(gca, 'FontSize', 12, 'FontName', 'Times', 'LineWidth', 1)

title('(a)', 'Units', 'normalized', 'Position', [0.04, 0.93, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_a.pdf'],'ContentType','vector');
end

figure;
bar(1:length(model), T.MeanPerStep, 'EdgeColor', 'none', 'FaceColor', '#6888F5');

hold on;
er = errorbar(T.MeanPerStep, T.StdPerStep, 'LineWidth', 1.5);
er.Color = [0 0 0];
er.LineStyle = 'none';

box off
grid on
xticklabels(model)
ylabel('Average simulation time of per step (s)')
set(gca, 'FontSize', 12, 'FontName', 'Times', 'LineWidth', 1)

title('(b)', 'Units', 'normalized', 'Position', [0.04, 0.93, 0], 'FontWeight', 'bold');

if isOutputToFile
    exportgraphics(gcf,[figpath, baseFileName,'_b.pdf'],'ContentType','vector');
end
