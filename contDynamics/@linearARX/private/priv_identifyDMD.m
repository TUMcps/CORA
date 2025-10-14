function sys = priv_identifyDMD(traj,p)
% priv_identifyDMD - Identifies a linear ARX system from trajectory data 
%   using Dynamic Mode Decomposition (DMD) according to [1]
%
% Syntax:
%    sys = priv_identifyDMD(traj,p)
%
% Inputs:
%    traj - trajectory data storing the states, times, and inputs for 
%             multiple trajectories
%    p - number of past time steps
%
% Outputs:
%    sys - identified linearARX object
%
% References:
%    [1] J.L. Proctor and et al. "Dynamic mode decomposition with control" 
%         SIAM Journal on Applied Dynamical Systems 15, 1 (2016), 142–161
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearARX/identify

% Authors:       Niklas Kochdumper
% Written:       30-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % split the data into single data points
    points = aux_getDataPoints(traj,p);

    % apply dynamic mode decomposition (DMD) for all ranks of the SVD
    Alist = aux_dynamicModeDecomposition(points.x,points.xNext,points.u);

    % select the matrix that best fits the data
    errBest = inf;
    dt = traj(1).t(2) - traj(1).t(1);
    n = size(points.xNext,1); m = size(traj(1).u,1);

    for i = 1:length(Alist)

        % construct linear ARX system object
        A = Alist{i}(:,1:n*p); B = Alist{i}(:,n*p+1:end);
        A_ = mat2cell(A,n,n*ones(1,p))';
        B_ = mat2cell(B,n,m*ones(1,p+1))';
        sysTemp = linearARX(A_,B_,dt);

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

function points = aux_getDataPoints(traj,p)
% transform the trajectory into a list of data points

    points.x = [];
    points.xNext = [];
    points.u = [];

    for i = 1:length(traj)

        for j = p+1:size(traj(i).x,2)
            points.xNext = [points.xNext traj(i).x(:,j)];
            points.x = [points.x, ...
                              reshape(fliplr(traj(i).x(:,j-p:j-1)),[],1)];

            if ~isempty(traj(i).u)
                points.u = [points.u, ...
                              reshape(fliplr(traj(i).u(:,j-p:j)),[],1)];
            end
        end
    end

    if isempty(points.x)
        throw(CORAerror("CORA:notSupported", ['Trajectories are too ', ... 
                                'short for identifying an ARX model!']));
    end
end

% ------------------------------ END OF CODE ------------------------------
