function sys = priv_identifySparse(traj,phi)
% priv_identifySparse - Identifies a nonlinear discrete-time system from 
%   trajectory data using the SINDy approach [1]
%
% Syntax:
%    sys = priv_identifySparse(traj,phi)
%
% Inputs:
%    traj - trajectory data 
%    phi - template functions phi(x,u) for the dynamics 
%          x(k+1) = A*phi(x(k),u(k)) (function handle)
%
% Outputs:
%    sys - identified nonlinearSysDT object
%
% References:
%    [1] S.L. Brunton and et al. "Discovering governing equations from data 
%         by sparse identification of nonlinear dynamical systems, 2016
%    [2] J.L. Proctor and et al. "Dynamic mode decomposition with control" 
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
    points = getDataPoints(traj, true, phi);

    % perform Lasso regression for each dimension
    A = zeros(size(points.xNext,1),size(points.phi,1));

    for i = 1:size(points.xNext,1)
        a = lasso(points.phi',points.xNext(i,:)');
        A(i,:) = a(:,1)';
    end

    % refine the sparse model using Dynamic Mode Decomposition [2]
    sys = aux_refineDynamicModeDecomposition(A,phi,traj,points);

    % construct nonlinear discrete-time system object
    %sys = nonlinearSysDT(@(x,u) A*phi(x,u),n,m);
end


% Auxiliary functions -----------------------------------------------------

function sys = aux_refineDynamicModeDecomposition(A,phi,traj,points)
% refine the sparse model using Dynamic Mode Decomposition (DMD)

    w = warning(); warning('off');

    % initialization
    n = size(traj(1).x,1);
    m = size(traj(1).u,1);

    dt = traj(1).t(2) - traj(1).t(1);

    % recompute optimal parameters for the model using Dynamic Mode
    % Decomposition [2]
    Anew = cell(size(A,1),1); max_rank = 0;

    for i = 1:size(A,1)
        phi_ = diag(abs(sign(A(i,:))))*points.phi;
        Anew{i} = aux_dynamicModeDecomposition(phi_,points.xNext(i,:));
        max_rank = max(max_rank,length(Anew{i}));
    end

    % loop over matrices with different rank parameter for the Singular
    % Value Decomposition and select rank that results in smallest error
    err = inf;

    x = sym('x',[n,1]); u = sym('u',[m,1]);
    phiOrig = phi(x,u);

    for i = 1:max_rank

        % construct coefficient matrix for this rank
        phi_ = [];
        for j = 1:length(Anew)
            ind = min(i,length(Anew{j}));
            phi_ = [phi_;Anew{j}{ind}*diag(abs(sign(A(j,:))))*phiOrig];
        end

        % construct nonlinear system object
        f = matlabFunction(phi_,'Vars',{x,u});

        sysTmp = nonlinearSysDT(f,dt,n,max(1,m));

        % compute error
        errTmp = computeError(traj,sysTmp);

        if errTmp < err
            err = errTmp;
            sys = sysTmp;
        end
    end

    warning(w);
end

function A = aux_dynamicModeDecomposition(X1,X2)
% compute the matrix X2 = A*X1 that best fits the data using the approach
% in Equation (2.7) in [2]

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
