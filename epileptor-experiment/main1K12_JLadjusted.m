resultsDir = setupEpileptorPaths();

tic

% Set the parameter tau0
tau0 = 6667; 

% Specify the patient identifier
patient = 'P1';

K_orig = loadPatientWeights(patient); 

% Normalize the connectivity matrix K
K = normal(K_orig);

% Get the size of the connectivity matrix (number of nodes N)
[N, ~] = size(K); 

% Initialize cell arrays to store results and node indices
results = cell(1, 85);
node_indices = cell(1, 85);

% Store the original connectivity matrix and node indices
results{1} = K;  
node_indices{1} = 1:N; 

% Loop to remove each node and construct modified connectivity matrices
%##### Are we supposed to normalise again?
for i = 1:N
    K_modified = K;  
    K_modified(i, :) = [];  % Remove the i-th row
    K_modified(:, i) = [];  % Remove the i-th column
    results{i+1} = K_modified;  % Store the modified matrix
    node_indices{i+1} = setdiff(1:N, i);  % Update node indices
end

% Set numerical display format to long
format long

% Initialize the x0 parameter range
x0_range1 = -2.3;

% Define the x0 parameter range for the epileptogenic zone (EZ) nodes
x0_EZ_range1 = -2.3:0.1:-1.5;
disp(x0_EZ_range1);

% Initialize the stability matrix with NaN values
StabilityMatrix = NaN(N, N); 

% Initialize an array to store EZ node values
EZ_values = [];

%#### Flipping loop order to allow us to get at x0 changes for no resection
% first.

% Specify the index of the connectivity matrix to use
for k_idx = 1:N+1 % ##### But index 1 is no resection, is this handled?
    % Get the current connectivity matrix
    K_current = results{k_idx};
    [N_current, ~] = size(K_current); 

    % Iterate over original EZ nodes
    for orig_EZ = 1:N
        fprintf('Processing original EZ = %d\n', orig_EZ);

        % Get the current node indices
        current_indices = node_indices{k_idx}; 
        
        % Check if the original EZ node is in the current node indices
        if ismember(orig_EZ, current_indices)
            % Find the position of the EZ node in the current indices
            EZ = find(current_indices == orig_EZ); 
            % Add the original EZ node to EZ_values
            EZ_values = [EZ_values, orig_EZ]; 
            % Initialize variables to record failed and last successful x0_EZ values
            failed_x0_EZ = NaN; 
            last_real_x0_EZ = x0_EZ_range1(1); %##### Was NaN, whcih is better?
            found_transition = false;

            % Iterate over the range of x0_EZ values
            for ix0_EZ = 1:length(x0_EZ_range1)
                % Initialize the x0 vector for all nodes
                x0 = x0_range1 + zeros(N_current, 1);
                % Set the x0 value for the EZ node
                x0(EZ) = x0_EZ_range1(ix0_EZ);
                % Initialize the Z0 vector
                Z0 = 3 + zeros(N_current, 1); 
                % Define the function handle for fsolve
                one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K_current, tau0);
                % Set options for fsolve
                opt = optimset('TolFun', 1e-14, 'TolX', 1e-14);
                % Use fsolve to solve for the fixed point
                [z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);

                fprintf('Exit flag: %d\n', exitflag);

                % Check if fsolve was successful
                if exitflag == 1
                    last_real_x0_EZ = x0_EZ_range1(ix0_EZ);
                else
                    warning('Failed to find fixed point, skipping current x0(EZ) value.');
                    failed_x0_EZ = x0_EZ_range1(ix0_EZ);
                    % Perform a finer search between the last successful and failed x0_EZ values
                    fine_x0_EZ_range = linspace(last_real_x0_EZ, failed_x0_EZ, 11); 
                    for fine_x0_EZ = fine_x0_EZ_range(2:end-1) 
                        x0(EZ) = fine_x0_EZ; 

                        % Recompute Z0 and solve for the fixed point
                        Z0 = 3 + zeros(N_current, 1);
                        one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K_current, tau0);
                        [z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);

                        % Check if fsolve was successful
                        if exitflag == 1
                            last_real_x0_EZ = fine_x0_EZ;
                            continue;
                        else
                            fprintf('Fine stability1 failed at x0(EZ) = %f\n', fine_x0_EZ);
                            failed_x0_EZ = fine_x0_EZ;
                            break;
                        end
                    end
                    %#### By definition, we must have found the transition now:
                    break;
                end
            end

            if (k_idx == 1)
                fprintf('\nEZ=%d, no resection, stable down to x0=%.3f\n\n', orig_EZ, failed_x0_EZ);
            else
                fprintf('\nEZ=%d, resected=%d, stable down to x0=%.3f\n\n', orig_EZ, k_idx-1, failed_x0_EZ);
            end
            % Record the failed x0_EZ value in the stability matrix
            StabilityMatrix(orig_EZ, k_idx) = failed_x0_EZ;
        end
    end

    if (k_idx == 1)
        % Save a matrix of the max x0 values for stability for no resection
        noResectionMaxX0ValuesForStability = StabilityMatrix(:, 1);
        save(fullfile(resultsDir, 'noResectionMaxX0ValuesForStability.mat'), 'noResectionMaxX0ValuesForStability');
    end
    
end

% Determine the range of processed EZ nodes
EZ_range_min = min(EZ_values);
EZ_range_max = max(EZ_values); 

% Generate the filename to save the stability matrix
filename = sprintf('StabilityMatrix_%d_to_%d.mat', EZ_range_min, EZ_range_max);

% Save the stability matrix to a .mat file
outFile = fullfile(resultsDir, filename);
save(outFile, 'StabilityMatrix');
disp(['Saved ', outFile, ' successfully.']);

toc
