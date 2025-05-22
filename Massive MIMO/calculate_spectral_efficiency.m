function C = calculate_spectral_efficiency(G, P, noise)
[p, ~, ~] = waterfilling(G, P, noise);
[~, S, ~] = svd(G, 'econ');
lambda = diag(S);
C_i = log2(1 + (p .* lambda.^2) / noise);
C = sum(C_i);
end