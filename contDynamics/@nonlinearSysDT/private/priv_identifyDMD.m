function [sys,A] = priv_identifyDMD(traj,phi)
% priv_identifyDMD - Identifies a nonlinear discrete-time system from 
%   trajectory data using Dynamic Mode Decomposition (DMD) according to [1]
%
% Syntax:
%    sys = priv_identifyDMD(traj,phi)
%
% Inputs:
%    traj - trajectory data
%    phi - template functions phi(x,u) for the dynamics 
%          x(k+1) = A*phi(x(k),u(k)) (function handle)
%
% Outputs:
%    sys - identified nonlinearSysDT object
%    A - identified parameter matrix
%
% References:
%    [1] J.L. Proctor and et al. "Dynamic mode decomposition with control" 
%         SIAM Journal on Applied Dynamical Systems 15, 1 (2016), 142–161
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT/identify

% Authors:       Niklas Kochdumper
% Written:       17-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % split the data into single data points
    points = getDataPoints(traj, false, phi);
    
    % apply dynamic mode decomposition (DMD) for all ranks of the SVD
    Alist = aux_dynamicModeDecomposition(points.phi,points.xNext);

    % select the matrix that best fits the data
    w = warning; warning('off');
    errBest = inf;
    n = size(traj(1).x,1); dt = traj(1).t(2) - traj(1).t(1);
    m = 1;
    if ~isempty(traj(1).u)
        m = size(traj(1).u,1);
    end
    for i = 1:length(Alist)

        % construct nonlinear system object
        sysTmp = nonlinearSysDT(@(x,u) Alist{i}*phi(x,u),dt,n,m);

        % compute the error between the system approximation and the real data
        err = computeError(traj,sysTmp);

        if err < errBest
            sys = sysTmp; A = Alist{i};
            errBest = err;
        end
    end

    warning(w);
end


% Auxiliary functions -----------------------------------------------------

function A = aux_dynamicModeDecomposition(X1,X2)
% compute the matrix X2 = A*X1 that best fits the data using the approach
% in Equation (2.7) in [1]

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
