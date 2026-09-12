function [x0_c, rho_c] = findCriticalX0(i, N, K, tau0, x0_base, x0_upper, x0_step, bisectTol, opt)
%FINDCRITICALX0  Smallest x0_i for node i at which rho(C) >= 1.
%
%   With all other nodes held at x0_j = x0_base, coarse forward sweep with
%   step x0_step until rho(C) >= 1 or fsolve fails (treated as seizure
%   onset). Then bisect between the last stable and first unstable x0_i
%   to tolerance bisectTol.

x0_vec_base = x0_base + zeros(N, 1);
Z0 = 3 + zeros(N, 1);

last_stable     = x0_base;
last_stable_rho = NaN;
last_stable_z   = Z0;
first_unstable  = NaN;
first_unstable_rho = NaN;

coarse_grid = x0_base:x0_step:x0_upper;
prev_z = Z0;
fsolve_failure_streak = 0;
for k = 1:length(coarse_grid)
    [stable, rho_val, z_out] = evalEpileptorStability(coarse_grid(k), i, x0_vec_base, ...
        K, tau0, prev_z, opt);
    if stable
        last_stable     = coarse_grid(k);
        last_stable_rho = rho_val;
        last_stable_z   = z_out;
        prev_z          = z_out;
        fsolve_failure_streak = 0;
    elseif ~isfinite(rho_val) && fsolve_failure_streak < 1
        % Tolerate a single isolated fsolve failure: skip without declaring
        % instability so a transient numerical glitch does not truncate the sweep.
        fsolve_failure_streak = fsolve_failure_streak + 1;
        continue;
    else
        first_unstable     = coarse_grid(k);
        first_unstable_rho = rho_val;
        break;
    end
end

if isnan(first_unstable)
    x0_c = NaN;
    rho_c = last_stable_rho;
    return;
end

lo = last_stable;
hi = first_unstable;
hi_rho = first_unstable_rho;
warm_z = last_stable_z;
while (hi - lo) > bisectTol
    mid = 0.5 * (lo + hi);
    [stable, rho_val, z_out] = evalEpileptorStability(mid, i, x0_vec_base, K, tau0, warm_z, opt);
    if stable
        lo = mid;
        warm_z = z_out;
    else
        hi = mid;
        hi_rho = rho_val;
    end
end

x0_c = hi;
rho_c = hi_rho;
end
