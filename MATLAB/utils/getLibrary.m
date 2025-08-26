function [Theta, Charset] = getLibrary(X, maxDegree)

len     = size(X, 1);        % data length
dim     = size(X, 2);        % data dimension
D       = maxDegree;         % max degree of terms


% % [1, 2, 3, 4, ... , n, n+1]个变量中最多选择k个
max_k   = D;
n       = dim;
% N       = (n+1)^max_k;      % 初始N值
% % N       = nchoosek(dim+D, D);
% mat     = zeros(max_k, N);
% num     = libpointer('double',1);
% k       = 1;
% temp    = zeros(max_k, 1);
% mat     = backTrack(n, N, mat, num, temp, k, max_k);
% for i = 1:N                 % 处理mat矩阵
%     mat(:, i) = sort(mat(:, i));
% end
% mat = unique(mat', 'rows')';
% N       = size(mat, 2);     % 更新N的值

mat = generate_polynomial_combinations(dim, D);
N       = size(mat, 1);     % 更新N的值



% d = max_polyorder;
% [~,dim] = size(X);    

% d = D;
% 
% % Generate all combined forms
% [Ind{1:dim}] = ndgrid(0:d);
% C = zeros((d+1)^dim, dim);
% for i = 1:dim
%     C(:, i) = Ind{i}(:);
% end
% 
% % Select the items that meet the criteria (order <=d) 
% ind_situable = sum(C, 2) <= d;
% C_situ = C(ind_situable, :);
% 
% % Order the items in a certain order
% mat = zeros(size(C_situ));
% ind_start = 1;
% for i = 0:d
%     order_i = sum(C_situ, 2)==i;
%     cnt_order_i = sum(order_i);
%     mat(ind_start:(ind_start+cnt_order_i-1), :) = C_situ(order_i, :);
%     ind_start = ind_start + cnt_order_i; 
% end
% N       = size(mat, 1); 
% % disp(mat)

var     = cell(1, n+1);     % 生成字符串形式变量名
for i = 1:n+1
   if i == 1
       var{i} = '1';
   else
       var{i} = ['x', num2str(i-1)]; 
   end
end

% disp(var)

Theta = ones(len, N);       % 生成观测矩阵/字典矩阵 Theta
% for i = 1:N
%    for j = 1:max_k
%        ind = mat(j, i);
%        if ind ~= 1
%            Theta(:, i) = Theta(:, i).*X(:,ind-1);
%        end
%    end
% end
for i = 1:N
    Theta(:, i) = prod(X.^mat(i,:),2);
end 

% Charset = cell(N, 1);
% for i = 1:N
%    if sum(mat(:, i)) == max_k
%       Charset{i} = '1'; 
%    else
%       var_str = '';
%       for j = 1:max_k
%          ind = mat(j, i);
%          if ind ~= 1
%             var_str = [var_str, var{ind}];
%          end 
%       end
%       Charset{i} = var_str;
%    end
% end
Charset = cell(N, 1);
for i = 1:N
   if sum(mat(i,:)) == 0
      Charset{i} = '1'; 
      % fprintf('Charset{%d}= %s \n', i, Charset{i});
   else
      var_str = '';
      for j = 1:dim
         ind = mat(i,j);
         % fprintf('i= %d \n', i);
         % if ind ~= 0
         %    var_str = [var_str, var{j+1}];
         % end
         for k = 1:ind
            var_str = [var_str, var{j+1}];
         end
         % if ind ~= 1
         %    var_str = [var_str, var{ind}];
         % end 
      end
      Charset{i} = var_str;
      % fprintf('Charset{%d}= %s \n', i, Charset{i});
   end
end
% disp(Charset);
% No = (1:length(Charset))';
% tbl = table(No, Charset);
% disp(tbl);



