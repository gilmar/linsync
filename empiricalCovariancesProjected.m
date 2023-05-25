function [empiricalCovarianceProjected] = empiricalCovariancesProjected(C, N, discreteTime, S, dt)
%
% Returns an empirical measurement of the covariance in the space orthogonal to the zero mode psi_0=[1 1 1 ... 1]
%  for either continuous time Ornstein-Uhelnbeck or discrete time AR
%  dynamics as required.
% 
% There was no original implementation of this in the Barnett ncomp
% toolkit, so we built the Cholesky decomposition here following their
% guidelines. Differently to that, we are projecting orthogonal to the
% fully synchronized mode at every time step, since we want the projected
% covariance matrix in the end anyway and this will give a more stable
% result.
%
% Inputs:
% - C - network coupling matrix which can be NxN for the no-delays case,
%      (where C(i,j) = i->j weight, in row vector form), or N(tau+1)xN(tau+1) for
%      the case of additional lagged connections up to tau (in which case
%      it is actually the B matrix
% - N - the size of the network for which we want to compute the covariance
%      of. This case be a subset of the supplied network C (for the lagged
%      case for instance)
% - discreteTime - whether the time-series is estimated from a discretized process (true) or use of exact method assuming continuous process (false).
% - S - number of time steps (sample size) of the time series to compute the covariance over.
% - dt - time interval for iterating the exact method (continuous time). Ignored for discrete time.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

NtimeTauPlus1_1 = size(C,1);
NtimeTauPlus1_2 = size(C,2);
if (NtimeTauPlus1_2 ~= NtimeTauPlus1_1) 
    error('C is not square');
end

% For system X:
I = eye(N);
% Averaging projection operator:
G = 1 ./ N * ones(N);
% Unaveraging operator:
U = I - G;

hasDelays = (N < NtimeTauPlus1_1);
if (hasDelays)
    % For embedded system Z:
    B = C; % For matching notation to text, we will now use B instead of C for the delays case
    zLength = size(B,1);
    tauPlus1 = zLength / N;
    IZ = eye(zLength);
    IX = zeros(zLength);
    IX(1:N,1:N) = I;
    GX = [repmat(G, 1, tauPlus1); zeros(N*(tauPlus1-1), N*tauPlus1)];
    UX = IZ - GX;
end

% Burn-in time before keeping samples for analysis:
initialTimeSteps = 1000 + S; % For no lambda_0=1, this parameter doesn't seem to matter much, it settles quickly. No difference using 1000 or 100000 steps. For having \lambda_0=1, having at least 10000 steps would be advisable, though the accuracy seems to depend more on S than this.

tol = 10; % Multiples of machine epsilon within which we want to consider a value to be zero.

if (discreteTime)
    % Run the discrete time multivariate AR (autoregressive) process

    if (~hasDelays)
        X = normrnd(0, 1, 1, N); % initialise row vector, with std deviation of the noise
        Xproj = X * U; % Projection away from the synchronized state
        allXproj = zeros(S,N); % Samples of X * U
        sampleNum = 1;
        for t = 1 : initialTimeSteps + S
            % Make the next time step:
            randNoise = normrnd(0, 1, 1, N);
            % Since we're going to project by U in the end anyway to compute sync, we're going
            % to make the projection at every step of the dynamics - this
            % avoids the values of X alone becoming unstable and compromising
            % the final empirical calculation.
            newXproj = Xproj * C * U + randNoise * U;
            if (t > initialTimeSteps)
                % Keep the new values
                allXproj(sampleNum, :) = newXproj;
                sampleNum = sampleNum + 1;
            end
            Xproj = newXproj;
        end
    else
        % Case for delays where X is embedded as system Z
        Z = normrnd(0, 1, 1, zLength); % initialise row vector, with std deviation of the noise
        Zproj = Z * UX; % Projection of Z away from the synchronized state of X
        allXproj = zeros(S,N); % Samples of X * UX
        sampleNum = 1;
        for t = 1 : initialTimeSteps + S
            % Make the next time step:
            randNoiseZ = [normrnd(0, 1, 1, N), zeros(1, zLength - N)];
            % Since we're going to project by UX in the end anyway to compute sync, we're going
            % to make the projection at every step of the dynamics - this
            % avoids the values of X alone becoming unstable and compromising
            % the final empirical calculation.
            newZproj = Zproj * B * UX + randNoiseZ * UX;
            if (t > initialTimeSteps)
                % Keep the new values for X (just the X part of Z)
                allXproj(sampleNum, :) = newZproj(1:N);
                sampleNum = sampleNum + 1;
            end
            Zproj = newZproj;
        end
    end
else
    % Run the Ornstein-Uhlenbeck process using exact method
    
    % How do the expected values of X adjust over dt: (appendix A of Barnett et al)
    %  <X(t+s)> = X(t) * meanMultiplier
    meanMultiplier = expm(-(I-C)*dt);

    % Generate the covariance matrix M(dt) for how the variables X covary
    %   after we make the dt time step if they started with no covariance:
    symmetricC = isempty(find(C - C' > tol*eps));
    [MsProjected,err] = contCon2LaggedCovProjected(C,dt,symmetricC,1000,tol,true);
    if (err ~= 0)
        return;
    end
    % Make a Cholesky decomposition - returns the upper triangle, or if MsProjected is rank deficient,
    %  which it is because we project orthogonal to the zero mode psi_0, then (from matlab help):
    % "R is an upper triangular matrix of size q-by-n so that the L-shaped region of the first q rows and first q
    %  columns of R'*R agree with those of A."
    [cholR,p] = chol(MsProjected);
    fprintf('chol returned p=%d, giving q=%d rows and cols in cholR\n', p, p - 1);
    if (p > 0)
        % Pad the cholesky decomposition with zeros to make an NxN matrix to multiply the noise by.
        % We assume that the zero multiplications here still have something come through for those columns
        %  because of the unaveraging procedure.
        noiseMultiplier = [cholR, zeros(p - 1, N - (p - 1)); zeros(N - (p - 1), N)];
    else
        noiseMultiplier = cholR;
    end

    X = normrnd(0, 1, 1, N) * sqrt(dt); % row vector like Barnett and Bullock, with std dev of the noise to start with
    Xproj = X * U;
    allXproj = zeros(S,N);
    sampleNum = 1;
    for t = 1 : initialTimeSteps + S
        % Make the next time step:
        % Since we're going to project by U in the end anyway to compute sync, we're going
        % to make the projection at every step of the dynamics - this
        % avoids the values of X alone becoming unstable and compromising
        % the final empirical calculation.

        % The mean of each vector after projection will be:
        newXProjMean = Xproj * meanMultiplier * U;
        % and the updates to the process over s will have covariance matrix MsProjected,
        %  so we work out the variances via a Cholesky decomposition (App A of Barnett):
        newXProj = newXProjMean + normrnd(0, 1, 1, N) * noiseMultiplier * U;
        if (t > initialTimeSteps)
            % Keep the new values
            % allX(:, t - multiplesOfSForSimulating*S) = newX; % If we're doing the S after 10*S
            allXproj(sampleNum, :) = newXProj;
            sampleNum = sampleNum + 1;
        end
        Xproj = newXProj;
    end
end

% If we had not made the projection at every step of the dynamics above, we would now as follows:
% empiricalCovarianceProjectedFromCov = cov(allX * U);
% But in theory the allX above may diverge (in practise we won't be using enough samples for
%  that to happen, but we'll get a better estimate this way anyway)
%  because the time series may have diverged before we can take the
%  covariance. So instead we project at every time step above, achieves
%  same as projecting now but is more stable.

empiricalCovarianceProjected = cov(allXproj);

end

