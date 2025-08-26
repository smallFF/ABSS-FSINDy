function file_name = generate_lorenz96(n, F)
file_name = ['lorenz96_n_',num2str(n),'_F_',num2str(F)];

fid = fopen([file_name,'.m'],'w');


fprintf(fid, 'function dX = %s(t, X)\n', file_name);
fprintf(fid, 'dX = zeros(%d, 1);\n', n);

% dx = (x_(i+1) - x_(i-2))*x_(i-1) - x_i + F;
% x_(i-n) = x_(i+n) = x_i;
for i = 1:n
    ind_i_plus_1 = modify_ind(i+1, n);
    ind_i_minus_2 = modify_ind(i-2, n);
    ind_i_minus_1 = modify_ind(i-1, n);
    dX_i = sprintf('dX(%d) = (X(%d) - X(%d))*X(%d) - X(%d) + %f;\n',...
        i, ind_i_plus_1, ind_i_minus_2, ind_i_minus_1, i,F);
    fprintf(fid, dX_i);
end
fclose(fid);

end


function ind = modify_ind(i, n)
    if i > n
        ind = i - n;
    elseif i < 1
        ind = i + n;
    else
        ind = i;
    end
end