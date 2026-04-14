%% set basic parameter
% D_st = (1/N)*trace(Omega) via covariancesGaussianNet; see README Section 1.6 and computeSyncResults.m
resultsDir = setupEpileptorPaths();
MaxK = 100000000;
tau0 = 6667;
patient = 'P1';
EZ = 8;
K = loadPatientWeights(patient);
K = normal(K); % Normalize
[N, ~] = size(K); % Get the number of nodes
%- create I_NEZ
I = eye(N);
I (EZ, EZ)=0;

%%- initialize
format long
x0_range1 = -2.30;
x0_EZ_range1 = -1.48:0.01:-1.47;
disp(x0_EZ_range1);
StabilityMatrix1 = NaN(length(x0_range1), length(x0_EZ_range1));
TotalCouplingMatrix = NaN(length(x0_range1), length(x0_EZ_range1), N, N);
LargestLambdaMatrix1 = NaN(length(x0_range1), length(x0_EZ_range1));
CILargestLambdaMatrix1 = NaN(length(x0_range1), length(x0_EZ_range1));
BeAwayStabilityMatrix = NaN(length(x0_range1), length(x0_EZ_range1), N);
PushAwayStabilityMatrix = NaN(length(x0_range1), length(x0_EZ_range1), N);
ExitFlagMatrix = NaN(length(x0_range1), length(x0_EZ_range1));
err1 = NaN(length(x0_range1), length(x0_EZ_range1));
err2 = NaN(length(x0_range1), length(x0_EZ_range1));

%%-dual loop x0 and x0(EZ)
for ix0 = 1:length(x0_range1)
   for ix0_EZ = 1:length(x0_EZ_range1)
       x0 = x0_range1(ix0) + zeros(N, 1);
       x0(EZ) = x0_EZ_range1(ix0_EZ);
      
       %%- solve fixed-point
       Z0 = 3 + zeros(N, 1);
       one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K, tau0);
       opt = optimset('TolFun', 1e-14, 'TolX', 1e-14);
       [z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);
       fprintf('Exit flag: %d\n', exitflag);
       ExitFlagMatrix(ix0, ix0_EZ) = exitflag;

	% Check exitflag. If less than or equal to 0, exit the loop.
	if exitflag <= 0
            fprintf('Exit flag <= 0 encountered. Stopping execution.\n, Its seizure mode!');
    	    break;
	end

       %%- calculate Coupling Matrix for numerical solution
       C1 = CouplingMatrix(z_fixed_numerical, K, tau0);
       TotalCouplingMatrix(ix0, ix0_EZ, :, :) = C1;
              

       %%- calculate Largest Eigenvalue (Lambda)
       lambda_max1 = max(eig(C1));
       LargestLambdaMatrix1(ix0, ix0_EZ) = lambda_max1;
       fprintf('Largest Lambda1: %.3f\n', lambda_max1);


       %%-  check for stationarity and convergence of the C
       if real(lambda_max1) > 1
           fprintf('Real part of largest eigenvalue of C > 1. Stopping execution.\n');
    	    break;
       end
	
       %if abs(lambda_max1) > 1
           %fprintf('Absolute value of largest eigenvalue of C > 1. Stopping execution.\n, Our math couldnt handle it.');
    	   %break;
       %end


       %%- D_st from full stationary covariance Omega (continuous-time OU); same as stability mode in computeSyncResults
       err = NaN;
       covariance1 = NaN(N);
       stability1 = NaN;
       try
           [~, covariance1, ~, err] = covariancesGaussianNet(C1, false, MaxK, false, 1);
           stability1 = trace(covariance1) / N;
       catch ME
           warning('runEZ1:covariancesGaussianNet', '%s', ME.message);
       end
       err1(ix0, ix0_EZ) = err;
     
       %%- check seizure type
       %if err == 2
    	   %CI_product = C1 * I;
           %eigenvalues_CI = eig(CI_product);
           %max_abs_eigenvalue = max(abs(eigenvalues_CI));
           %if max_abs_eigenvalue > 1
               %disp('Its spreading seizure');
           %else
               %disp('Its non-spreading seizure');
           %end
       %end
	
       CI_product = C1 * I;
       eigenvalues_CI = eig(CI_product);
       max_abs_eigenvalue = max(abs(eigenvalues_CI));
       CILargestLambdaMatrix1(ix0, ix0_EZ) = max_abs_eigenvalue ;

       if max_abs_eigenvalue > 1
          disp('Its spreading seizure');
       else
          disp('Its non-spreading seizure');
       end

       %%- save D_st (stored in StabilityMatrix1 for backward compatibility)
       StabilityMatrix1(ix0, ix0_EZ) = stability1;
       fprintf('D_st (stability): %.3f\n', stability1);
      
       %%- per-node diagonal of Omega (C)
       BeAwayStability = NaN(1, N);
       omega_diag = diag(covariance1);
       for i = 1:N
           BeAwayStability(i) = omega_diag(i);
       end
       BeAwayStabilityMatrix(ix0, ix0_EZ, :) = BeAwayStability;
       fprintf('BeAwayStability: %.3f\n', BeAwayStability(EZ));

       %%- per-node diagonal of Omega for dynamics with coupling C'
       err2(ix0, ix0_EZ) = NaN;
       try
           [~, covariance_transpose, ~, err_t] = covariancesGaussianNet(C1', false, MaxK, false, 1);
           err2(ix0, ix0_EZ) = err_t;
           PushAwayStability = NaN(1, N);
           omega_t_diag = diag(covariance_transpose);
           for i = 1:N
               PushAwayStability(i) = omega_t_diag(i);
           end
           PushAwayStabilityMatrix(ix0, ix0_EZ, :) = PushAwayStability;
       catch ME
           warning('runEZ1:covariancesGaussianNetTranspose', '%s', ME.message);
       end
   end
end

%% Store all matrices in a single structure
results.StabilityMatrix1 = StabilityMatrix1;
results.TotalCouplingMatrix = TotalCouplingMatrix;
results.LargestLambdaMatrix1 = LargestLambdaMatrix1;
results.BeAwayStabilityMatrix = BeAwayStabilityMatrix;
results.PushAwayStabilityMatrix = PushAwayStabilityMatrix;
results.ExitFlagMatrix = ExitFlagMatrix;
results.err1 = err1;
results.err2 = err2;

%% Create a filename based on EZ number
filename = sprintf('EZ%d_results.mat', EZ);

%% Save the structure to a .mat file with the EZ-based name
save(fullfile(resultsDir, filename), 'results');

