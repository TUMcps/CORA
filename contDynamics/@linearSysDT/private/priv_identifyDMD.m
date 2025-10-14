function sys = priv_identifyDMD(traj)
% priv_identifyDMD - Identifies a linear discrete-time system from 
%   trajectory data using Dynamic Mode Decomposition (DMD) according to [1]
%
% Syntax:
%    sys = priv_identifyDMD(traj)
%
% Inputs:
%    traj - trajectory data storing the states, times, and inputs for 
%             multiple trajectories
%
% Outputs:
%    sys - identified linearSysDT object
%
% References:
%    [1] J.L. Proctor and et al. "Dynamic mode decomposition with control" 
%         SIAM Journal on Applied Dynamical Systems 15, 1 (2016), 142–161
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT/identify

% Authors:       Niklas Kochdumper
% Written:       17-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % split the data into single data points
    points = getDataPoints(traj, false);

    % apply dynamic mode decomposition (DMD) for all ranks of the SVD
    Alist = aux_dynamicModeDecomposition(points.x,points.xNext,points.u);

    % select the matrix that best fits the data
    errBest = inf;
    dt = traj(1).t(2) - traj(1).t(1);

    for i = 1:length(Alist)

        % construct linear system object
        n = size(Alist{i},1);
        A = Alist{i}(:,1:n); c = Alist{i}(:,n+1); B = Alist{i}(:,n+2:end);
        sysTemp = linearSysDT(A,B,c,dt);

        % compute the error between the system approximation and the real data
        err = computeError(traj,sysTemp);
        if err < errBest
            sys = sysTemp;
            errBest = err;
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function A = aux_dynamicModeDecomposition(X1,X2,U)
% compute the matrix X2 = A*X1 that best fits the data using the approach
% in Equation (2.7) in [1]

    X1 = [X1; ones(1,size(X1,2))];

    if ~isempty(U)
        X1 = [X1; U];
    end

    % singular value decomposition
    [V,S,W] = svd(X1,'econ');
    
    % construct matrices with different rank
    rankMax = sum(diag(S) > 0);
    A = cell(rankMax,1);

    for rank = 1:rankMax

        % reduce rank by removing the smallest singular values
        if ~isempty(rank) && rank < size(S,1)
            V_ = V(:,1:rank); S_ = S(1:rank,1:rank); W_ = W(:,1:rank);
        else
            V_ = V; S_ = S; W_ = W;
        end
    
        % compute resulting system matrix A
        A{rank} = X2*W_*diag(1./diag(S_))*V_';
    end
end

% ------------------------------ END OF CODE ------------------------------
