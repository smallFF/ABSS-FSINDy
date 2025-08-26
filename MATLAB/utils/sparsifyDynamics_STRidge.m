function Xi = sparsifyDynamics_STRidge(Theta,dXdt, threshold, rho)
if nargin < 4
    rho = 1e-6;
end

n = size(dXdt, 2);
A = Theta;
AtA = Theta'*Theta;
I = eye(size(AtA));

Xi = (AtA + rho*I) \ (A'*dXdt);

k_max = 20;

for k=1:k_max
    smallinds = (abs(Xi)<threshold);   % find small coefficients
    Xi(smallinds)=0;                % and threshold
    for ind = 1:n                   % n is state dimension
        biginds = ~smallinds(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        A_tmp = Theta(:, biginds);
        AtA_tmp = A_tmp'*A_tmp;
        I = eye(size(AtA_tmp));

        Xi(biginds,ind) = (AtA_tmp + rho*I) \ (A_tmp'*dXdt(:,ind));
    end
end