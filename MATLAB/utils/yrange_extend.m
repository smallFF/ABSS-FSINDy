function yrange_new= yrange_extend(y, p1, p2)
if nargin < 2
    % p1 = 0.25;
    % p2 = 0.25;
    p1 = 0.75;
    p2 = 0.5;
elseif nargin < 3
    p2 = 0.5; 
end

% y_max = max(y);
% y_min = min(y);
% y_max = y_max + abs(y_max)*p1;
% y_min = y_min - abs(y_min)*p2;

r_max = max( y );
r_min = min( y );
% if abs(r_max) > 0.5
%     r_max = r_max*(1+p1);
% end
% if abs(r_max - r_min) < 0.5
%     r_min = r_min - abs(r_max - r_min)*p2;
%     r_max = r_max + abs(r_max - r_min)*p2;
% elseif abs(r_max - r_min) > 0.5
%     r_min = - 1;
% end
margin = r_max - r_min;
r_max = r_max + margin*p1;
r_min = r_min - margin*p2;
yrange_new = [r_min, r_max];