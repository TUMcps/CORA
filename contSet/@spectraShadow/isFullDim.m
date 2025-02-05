function [res,subspace] = isFullDim(SpS)
% isFullDim - checks if a spectrahedral shadow is full-dimensional
%
% Syntax:
%    res = isFullDim(SpS)
%    [res,subspace] = isFullDim(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    res - true/false
%    subspace - (optional) Returns a set of orthogonal unit vectors
%               x_1,...,x_k such that S is strictly contained in
%               p+span(x_1,...,x_k)
%               where p is some arbitrary point in S
%               (here, 'strictly' means that k is minimal).
%               Note that if S is just a point, subspace=[].
%
% Example:
%    SpS1 = spectraShadow([eye(2) eye(2)]);
%    SpS2 = spectraShadow([1 0 -1 0;0 -1 0 1]);
%    isFullDim(SpS1) % true
%    isFullDim(SpS2) % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope/isFullDim

% Authors:       Adrian Kulmburg
% Written:       12-August-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimension
n = dim(SpS);

% If the answer has already been computed, just give out the answer
if ~isempty(SpS.fullDim.val)
    res = SpS.fullDim.val;
    if res
        subspace = eye(n);
        return
    else
        % Can't say anything about the subspace
    end
end

% check if the spectrahedral shadow is completely empty
if representsa_(SpS,'emptySet',1e-8)
    res = false;
    subspace = [];
    
    SpS.fullDim.val = false;
    SpS.bounded.val = true;
    SpS.emptySet.val = true;
    return
end

% Another special case: G is the zero matrix
if norm(SpS.G, 'fro') < 1e-8
    res = false;
    subspace = [];

    SpS.fullDim.val = false;
    SpS.bounded.val = true;
    SpS.emptySet.val = true;
    return
end
    
% Let's first check whether the spectrahedral shadow is empty, by searching
% for one point that is contained in P:
x0 = priv_findFeasiblePointSpectrahedron(SpS);
if isempty(x0)
    res = false;
    subspace = [];
    
    SpS.fullDim.val = false;
    SpS.bounded.val = true;
    SpS.emptySet.val = true;
    return
end

x0 = SpS.c + SpS.G * x0;

% Now, assuming x0 exists, we translate S by -x0 so that we can
% guarantee that the resulting polytope contains the origin:
SpS_iter = SpS - x0;

% Setup the list of vectors we seek
subspace = [];

for i = 1:dim(SpS)
    % Search for vectors in S_iter that are perpendicular to the ones
    % we have found so far
    x_iter = aux_maxNormPerpendicularSpectrahedron(SpS_iter,subspace);
    % If the solution is x_iter=0, then this means that there are no
    % non-trivial vectors that are perpendicular to those we have
    % found, and so we can stop the iteration
    if norm(x_iter,Inf) <= 1e-6
        break;
    end

    % test if optimize result is within tolerance for 0
    % to avoid numerical issues
    for j = 1:length(x_iter)
        if withinTol(x_iter(j),0,1e-8)
            x_iter(j) = 0;
        end
    end

    x_iter_unit = x_iter/norm(x_iter);

    % Now, it suffices to add x_iter to the list
    subspace = [subspace,x_iter_unit];

end

% Now that we have constructed the subspace we need, it is time to
% check what dimension it has; only if it is full dimensional, is the
% spectrahedral shadow non-degenerate
k = size(subspace,2);
if k == dim(SpS)
    res = true;
    % If that is the case, we can assume that subspace is the entire
    % space, and so a much easier ONB would be the following
    subspace = eye(dim(SpS));
else
    res = false;
end
SpS.fullDim.val = res;
    
end


% Auxiliary functions -----------------------------------------------------

function x = aux_maxNormPerpendicularSpectrahedron(SpS,X)
    % For a list of vectors X = [x_1,...,x_k], solve the problem
    %   max ||x||_oo
    %   s.t. x + Xw \in S,
    %        forall i: x_i'*x = 0
    %
    % This is equivalent to
    %   max_y max_x y'*x,
    %   s.t. x + Xw \in S,
    %        forall i: x_i'*x = 0,
    %        and where y is iterated over all standard vectors +-e_i.
    
    % read out dimension of the spectrahedral shadow
    n = dim(SpS);
    
    G = SpS.G;
    c = SpS.c;
    [A0,Ai] = priv_getCoeffMatrices(SpS);
    
    % Special case: If the Ai matrices are all zero AND X is empty, Yalmip
    % will have trouble for some reason to define the constraints, so we
    % need to manually compute the result here
    if isempty(X)
        all_Ai_zero = true;
        for i = 1:size(G,2)
            if ~all(all(~Ai{i}))
                all_Ai_zero = false;
                break
            end
        end
        if all_Ai_zero
            % The output now depends on whether or not A0 is PSD...
            if all(eig(A0)>=0)
                % If yes, then set of solutions is unbounded, and so we can
                % take any vector
                x = ones([dim(SpS) 1]);
            else
                % If not, the problem is infeasible
                x = [];
            end
            return
        end
    end
    
    % store maximum objective value and maximizer
    maximum = 0;
    maximizer = [];
    
    persistent options
    if isempty(options)
        if isSolverInstalled('mosek')
            options = sdpsettings('solver','mosek','verbose',0,'allownonconvex',0,'cachesolvers',1);
        else
            options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0,'cachesolvers',1);
        end
    end

    % loop over all dimensions (for y)
    for i=1:n
        beta = sdpvar(size(G,2),1);
        w = sdpvar(size(X,2),1);
        constraints = [];
        temp = A0;
        for j=1:size(G,2)
            temp = temp + beta(j) * Ai{j};
        end
        constraints = [constraints temp >= 0];

        if isempty(X)
            x = G*beta + c;
        else
            x = G*beta + c - X*w;
        end

        for j=1:size(X,2)
            constraints = [constraints X(:,j)'*x==0];
        end

        % We also limit the range of all x(i) to 1, to deal with the case
        % where the spectrahedron is unbounded
        % (Note that checking boundedness is expensive, so it's cheaper
        % overall to just limit the spectrahedron right out of the gate)
        for ell=1:size(G,1)
            constraints = [constraints x(ell)<=1 -x(ell)<=1];
        end
        
        % compute maximum for y = +e_i        
        cost = -x(i);
        
        yalmipDiagnostics = optimize(constraints,cost,options);

        
        val = value(cost);
        
        if val > maximum
            maximum = val;
            maximizer = value(x);
        end
        
        % Do the same for y = -e_i
        
        % compute maximum for y = -e_i        
        cost = x(i);
        
        yalmipDiagnostics = optimize(constraints,cost,options);
        
        val = -value(cost);
        
        
        if val > maximum
            maximum = val;
            maximizer = value(x);            
        end
    end

    % return maximizer
    x = maximizer;
end

% ------------------------------ END OF CODE ------------------------------
