Ps = zeros(S, 1);
lam_sq = diag(lam).^2;
inverse_gains = noise ./ lam_sq;
finite_inverse_gains = inverse_gains(isfinite(inverse_gains));

if isempty(finite_inverse_gains)
    disp('All channel gains are zero, and power cannot be allocated.');
    Ps = zeros(size(inverse_gains));
    return;
end

tau_low = min(finite_inverse_gains);
tau_high = max(finite_inverse_gains) + P;
tolerance = 1e-6;
max_iterations = 100;
current_iteration = 0;
optimal_tau = 0;

while (tau_high - tau_low) > tolerance && current_iteration < max_iterations
    tau_mid = (tau_low + tau_high) / 2;
    current_allocated_ps_temp = max(0, tau_mid - inverse_gains);
    current_P_sum = sum(current_allocated_ps_temp);
    
    if abs(current_P_sum - P) < tolerance
        optimal_tau = tau_mid;
        break;
    elseif current_P_sum < P 
        tau_low = tau_mid;
    else
        tau_high = tau_mid;
    end
    
    current_iteration = current_iteration + 1;
end

if optimal_tau == 0
    optimal_tau = (tau_low + tau_high) / 2;
end

Ps = max(0, optimal_tau - inverse_gains);

Ps = Ps * min(1, P / sum(Ps));

disp('--- Result of water injection method power distribution ---');
disp('Allocate power:');
disp(Ps);
disp(['The total sum of distributed power: ', num2str(sum(Ps))]);