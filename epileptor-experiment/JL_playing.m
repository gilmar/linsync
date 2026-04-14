% JL playing:
% D_st = (1/N)*trace(Omega) via covariancesGaussianNet; see README Section 1.6 and computeSyncResults.m
resultsDir = setupEpileptorPaths();
tic

MaxK = 10000000;

% Set the parameter tau0
tau0 = 6667; 

% Specify the patient identifier
patient = 'P1';

K_orig = loadPatientWeights(patient); 

%#### What does this do? In the code from TVB people
% Normalize the connectivity matrix K
K = normal(K_orig);

% Get the size of the connectivity matrix (number of nodes N)
[N, ~] = size(K); 

figure(1);
imagesc(K);
colorbar;
title('K matrix');

% Initialize the x0 parameter range
x0 = -2.3 + zeros(N, 1);
% Set up an EZ?:
% EZ = 64;
% x0(EZ) = -1.8;
% Initialize the Z0 vector
Z0 = 3 + zeros(N, 1); 

%#### 1D epileptor should use matrix equations to be faster.
% Define the function handle for fsolve
one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K, tau0);
% Set options for fsolve
opt = optimset('TolFun', 1e-14, 'TolX', 1e-14);
% Use fsolve to solve for the fixed point
[z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);

fprintf('Exit flag: %d\n', exitflag);

% Check if fsolve was successful
%######### Need to change the check in main code to match this I think:
if exitflag ~= 1
    warning('Failed to find fixed point. (error code %d)', exitflag);
    return;
end

%#### Check use of tau0 here:
% Compute the coupling matrix C1
C1 = CouplingMatrix(z_fixed_numerical, K, tau0);

figure(2);
imagesc(real(C1)); %% Does it matter if C1 is not real??
colorbar;
title('C matrix');

err = NaN;
covariance1 = NaN(N);
stability1 = NaN;
try
    [~, covariance1, ~, err] = covariancesGaussianNet(C1, false, MaxK, false, 1);
    stability1 = trace(covariance1) / N;
catch ME
    warning('Error computing D_st: %s', ME.message);
end

if (~isnan(err) && err ~= 0)
    warning('covariancesGaussianNet / con2cov reported err = %d\n', err);
end

fprintf('D_st (deviation from stability): %.4f\n', stability1);

figure(3);
sigma_i = diag(covariance1);
plot(sigma_i, '-x')
xlabel('Node index');
ylabel('\sigma_i^2');
title('Deviation from stability');

% Omega for coupling C' (transpose dynamics); compare D_st to stability1 above
err_t = NaN;
covariance1Transpose = NaN(N);
stability1Transpose = NaN;
try
    [~, covariance1Transpose, ~, err_t] = covariancesGaussianNet(C1', false, MaxK, false, 1);
    stability1Transpose = trace(covariance1Transpose) / N;
catch ME
    warning('Error computing D_st for C'': %s', ME.message);
end
if (~isnan(err_t) && err_t ~= 0)
    warning('covariancesGaussianNet (C'') reported err = %d\n', err_t);
end

fprintf('D_st for C transpose: %.4f\n', stability1Transpose);

figure(4);
sigma_i_driving = diag(covariance1Transpose);
plot(sigma_i_driving, '-x')
xlabel('Node index');
ylabel('\sigma_i^2');
title('Driving contribution to deviation from stability ');

%% Now do virtual resections on C matrix:
deltaStabilityResectionsC = zeros(1,N);
for resectedNode = 1 : N
    C_modified = C1;
    C_modified(resectedNode, :) = [];  % Remove the i-th row
    C_modified(:, resectedNode) = [];  % Remove the i-th column    
    err = NaN;
    try
        [~, covarianceResected, ~, err] = covariancesGaussianNet(C_modified, false, MaxK, false, 1);
        Nres = size(covarianceResected, 1);
        stabilityResected = trace(covarianceResected) / Nres;
    catch ME
        warning('Error computing D_st after C resection: %s', ME.message);
        stabilityResected = NaN;
    end
    if (~isnan(err) && err ~= 0)
        warning('covariancesGaussianNet (resected C) err = %d\n', err);
    end
    fprintf('Deviation from stability with %d resected is: %.4f\n', resectedNode, stabilityResected);
    % Delta: positive means it got more stable (deviated less after
    % resection)
    deltaStabilityResectionsC(resectedNode) = stability1 - stabilityResected;
end
figure(5);
plot(deltaStabilityResectionsC, '-x')
xlabel('Node index');
ylabel('\delta sigma_i^2');
title('Delta deviation from stability after resection C');

%% Now do virtual resections on K matrix, and normalise after? Is this the way?
deltaStabilityResectionsKNormed = zeros(1,N);
for resectedNode = 1 : N
    K_modified = K_orig;
    K_modified(resectedNode, :) = [];  % Remove the i-th row
    K_modified(:, resectedNode) = [];  % Remove the i-th column    
    % Normalize the connectivity matrix K_modified
    K_modified_normed = normal(K_modified);
    % Initialize the x0 parameter range, without an EZ
    x0 = -2.3 + zeros(N-1, 1);
    % Initialize the Z0 vector
    Z0 = 3 + zeros(N-1, 1); 
    % Define the function handle for fsolve
    one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K_modified_normed, tau0);
    % Set options for fsolve
    opt = optimset('TolFun', 1e-14, 'TolX', 1e-14);
    % Use fsolve to solve for the fixed point
    [z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);
    % fprintf('Exit flag: %d\n', exitflag);
    % Check if fsolve was successful
    % %% Need to change the check here:
    if exitflag ~= 1
        warning('Failed to find fixed point. (error code %d)', exitflag);
        return;
    end
    % Compute the coupling matrix C1
    C1_resected = CouplingMatrix(z_fixed_numerical, K_modified_normed, tau0);
    err = NaN;
    try
        [~, covarianceResected, ~, err] = covariancesGaussianNet(C1_resected, false, MaxK, false, 1);
        Nres = size(covarianceResected, 1);
        stabilityResected = trace(covarianceResected) / Nres;
    catch ME
        warning('Error computing D_st after K resection: %s', ME.message);
        stabilityResected = NaN;
    end
    if (~isnan(err) && err ~= 0)
        warning('covariancesGaussianNet (resected K) err = %d\n', err);
    end
    fprintf('Deviation from stability with %d resected is: %.4f\n', resectedNode, stabilityResected);
    % Delta: positive means it got more stable (deviated less after
    % resection)
    deltaStabilityResectionsKNormed(resectedNode) = stability1 - stabilityResected;
end

figure(6);
plot(deltaStabilityResectionsKNormed, '-x')
xlabel('Node index');
ylabel('\delta sigma_i^2');
title('Delta deviation from stability after resection K normed');

%% Save the results:
save(fullfile(resultsDir, 'JL_playing_results.mat'), 'deltaStabilityResectionsKNormed', 'deltaStabilityResectionsC', 'K', 'C1', 'sigma_i', 'sigma_i_driving');

toc
