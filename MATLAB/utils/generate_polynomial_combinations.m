function combinations = generate_polynomial_combinations(n, d)
    combinations = zeros(0, n); 

    for k = 0:d
        if k == 0
            combinations = [combinations; zeros(1, n)];
        else
            if n == 1
                combinations = [combinations; k];
            else
                total_pos = k + n - 1;

                s = nchoosek(1:total_pos, n-1);
                num_combs = size(s, 1);
                combs = zeros(num_combs, n);
                
                for i = 1:num_combs
                    current_s = s(i, :);

                    ext_s = [0, current_s, total_pos + 1];

                    diffs = diff(ext_s);
                    
                    a = diffs - 1;
                    combs(i, :) = fliplr(a);
                end
                combinations = [combinations; combs];
            end
        end
    end
end