function [A, rho] = correlation_analysis(Theta, Charset, type)
if nargin < 3
    type = 1;
end
if type == 1
    value = 'Pearson';
elseif type == 2
    value = 'Kendall';
elseif type == 3
    value = 'Spearman';
else
    error('Supportted type value: 1, 2, 3.');
end
rho = corr(Theta(:, 2:end), 'Type', value);
string_name = Charset(2:end);
for i = 1:length(string_name)
    s = string_name{i};
    string_name{i} = regexprep(s, '\d+', '_$0');
end

etea = 0.9;
D = size(rho, 1);
figure('Units', 'centimeters', 'Position', [3 3 17 15]);
h = heatmap(string_name, string_name, rho);
set(gca, 'FontSize', 15);

rho=abs(rho);
rho_1=rho.*triu(ones(D,D), 1);
value=find(rho_1>etea);
Num=length(value);

A=zeros(Num-D, 3);
count = 0;
for i=1:D-1
    for j=i+1:D
        if rho_1(i, j) > etea
            count = count + 1;
            A(count,:)=[i+1, j+1, rho_1(i, j)];

            fprintf(['%2d |Two dimensions with strong linear correlation are: The %dth dimension: ...' ...
                '%s and the %dth dimension: %s, and their correlation coefficient is:%f\n'], ...
                count, i+1, string_name{i}, j+1, string_name{j}, rho_1(i, j));
        end

    end
end

