function t= FPT(T, s1, s2, epsilon)
t = T(end);

% gamma = 0;
% ind = abs(s1-s2)./(gamma+abs(s1))>=epsilon;
ind = abs(s1-s2) >= epsilon;
t = T(find(ind, 1));
% if ind
%     t = T(find(ind, 1));
% else
%     t = T(end);
% end
