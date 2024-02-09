function [res,subspace] = isFullDim(P)
% isFullDim - checks if a polytope is full-dimensional
%
% Syntax:
%    res = isFullDim(P)
%    [res,subspace] = isFullDim(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    res - true/false
%    subspace - (optional) Returns a set of orthogonal unit vectors
%               x_1,...,x_k such that P is strictly contained in
%               center(P)+span(x_1,...,x_k)
%               (here, 'strictly' means that k is minimal).
%               Note that if P is just a point, subspace=[].
%
% Example:
%    P1 = polytope([1 0;0 1;-1 0; 0 -1],[1;1;0;0]);
%    P2 = polytope([1 0;0 1;-1 0; 0 -1],[1;0;0;0]);
%    isFullDim(P1)
%    isFullDim(P2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Authors:       Niklas Kochdumper, Viktor Kotsev, Mark Wetzlinger, Adrian Kulmburg
% Written:       02-January-2020 
% Last update:   08-August-2022
%                05-December-2022 (VK, new method)
%                14-December-2022 (MW, add call to MOSEK)
%                25-May-2023 (AK, added method to compute subspace)
%                28-May-2023 (MW, add quick check for pairwise annihilation)
%                27-July-2023 (MW, add 1D method)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% dimension
n = dim(P);

% check if fullDim property is already set
if ~isempty(P.fullDim.val)
    res = P.fullDim.val;
    if nargout == 1
        return
    else
        % continue and compute subspace below...
    end
end

% if polytope has V-repesentation, then check rank of vertex matrix
if ~isempty(P.V.val)
    rankV = rank(P.V.val,1e-10);
    res = rankV == n;
    if res
        subspace = eye(n);
    else
        [Q,R] = qr(P.V.val);
        subspace = Q(:,rankV);
    end
    P.fullDim.val = res;
    return
end

% check if the polytope is completely empty
if representsa_(P,'emptySet',1e-10)
    res = false;
    subspace = [];
    return
end

% simpler version for 1D polytopes
if n == 1
    [res,subspace] = aux_1D(P);
    
    % set property
    P.fullDim.val = res;
    return
end


if nargout < 2
    % If the user only wants to know whether the polytope is degenerate, we
    % can check this rapidly using linear programming

    % Existence of equality constraints implies degenerate case (set may
    % also be empty)
    if ~isempty(P.Ae)
        res = false;
        P.fullDim.val = false;
        return
    end

    % number of inequality constraints and dimension
    [nrCon,n] = size(P.A);

    % normalize constraints
    P_norm = normalizeConstraints(P,'A');

    % quick check: any pair of inequalities of type
    %   ax <= b, ax >= b_ (-ax <= -b_), with b_ >= b?
    dotprod = tril(P_norm.A * P_norm.A');
    idxAntiParallel = find(withinTol(dotprod,-1));
    for i=1:length(idxAntiParallel)
        % read out matrix indexing from linear indexing
        [r,c] = ind2sub([nrCon,nrCon],idxAntiParallel(i));
        % check offset values, two perspectives:
        %  ax <= b1, -ax <= b2 <=> ax >= -b2  ->  -b2 >= b1?
        % -ax <= b1 <=> ax >= -b1, ax <= b2   ->  -b1 >= b2 <=> -b2 >= b1?
        % ...the one with the reversed sign has to be larger (or equal)!
        if -P_norm.b(c) > P_norm.b(r) || withinTol(-P_norm.b(c),P_norm.b(r))
            res = false;
            P.fullDim.val = false;
            return
        end
    end
    

    % solver info
    persistent isMosek
    if isempty(isMosek)
        isMosek = isSolverInstalled('mosek');
    end

    % 2-Norm of each row (always 1 since we normalized the constraints)
    A_norm = ones(nrCon,1);

    % extend inequality and equality constraints by one column
    A = [P.A A_norm];

    % cost function for linear program: minimize 2-norm of constraints
    f = [zeros(n,1); -1];

    % different solvers
    if isMosek

        % call MOSEK
        resMOSEK = msklpopt(f,A,-Inf(nrCon,1),P.b,...
            [-Inf(n,1);0],[],[],'minimize echo(0)');

        % read out solution
        if strcmp(resMOSEK.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
            % check if point is strictly interior (radius of center = 0)
            if withinTol(0,resMOSEK.sol.itr.dobjval,1e-9)
                res = false;
                P.fullDim.val = false;
            else
                res = true;
                P.fullDim.val = true;
            end
        elseif strcmp(resMOSEK.sol.itr.prosta,'PRIMAL_INFEASIBLE')
            % set is empty
            res = false;
            P.fullDim.val = false;
            P.emptySet.val = true;
        elseif strcmp(resMOSEK.sol.itr.prosta,'DUAL_INFEASIBLE')
            % numerical issues seem to come up only in the case of
            % full-dimensional unbounded polytopes
            res = true;
            P.fullDim.val = true;
            %throw(CORAerror('CORA:solverIssue'));
        end

    else
        % MATLAB linprog

        % linear program options
        persistent options
        if isempty(options)
            options = optimoptions('linprog','display','off');
        end

        % Solve Linear Program
        [~,r,exitflag] = linprog(f,A,P.b,[],[],[-Inf(n,1);0],[],options);

        if exitflag == 1
            % check if point is strictly interior
            if withinTol(0,r,1e-9)
                % r is 0 -> no ball fits in polytope -> degenerate
                res = false;
                P.fullDim.val = false;
            else
                % -r < 0 => radius of ball greater than 0 that fits in
                % polytope -> non-degenerate
                res = true;
                P.fullDim.val = true;
            end
        elseif exitflag == -2
            % set is empty
            res = false;
            P.fullDim.val = false;
            P.emptySet.val = true;
        elseif exitflag == -3
            % numerical issues seem to come up only in the case of
            % full-dimensional unbounded polytopes
            res = true;
            P.fullDim.val = true;
        elseif exitflag < 0
            throw(CORAerror('CORA:solverIssue'));
        end

    end

elseif nargout == 2
    % If on the other hand, the user is interested in the subspace, we have
    % to do the computations from above in a more specific way.
    
    % Let's first check whether the polytope is empty, by searching for one
    % point that is contained in P:
    x0 = aux_maxNormPerpendicularPolytope(P,[]);
    if isempty(x0)
        res = false;
        P.fullDim.val = false;
        P.emptySet.val = true;
        subspace = [];
        return
    end
    
    % Now, assuming x0 exists, we translate P by -x0 so that we can
    % guarantee that the resulting polytope contains the origin:
    P_iter = P - x0;
    
    % Setup the list of vectors we seek
    subspace = [];
    
    for i = 1:dim(P)
        % Search for vectors in P_iter that are perpendicular to the ones
        % we have found so far
        x_iter = aux_maxNormPerpendicularPolytope(P_iter,subspace);
        % If the solution is x_iter=0 (within same tolerance as used for
        % solving LPs), then this means that there are no non-trivial
        % vectors that are perpendicular to those we have found, and so we
        % can stop the iteration
        if norm(x_iter,Inf) <= 1e-8
            break;
        end
        
        % If we did find a new vector, we need to choose another vector,
        % along x_iter, that does not lie on a vertex. To do so, we will
        % also need to re-center P_iter
        x_iter_unit = x_iter ./ vecnorm(x_iter);
        % First, we compute how far we can scale x_iter with respect to the
        % origin, in the direction x_iter and -x_iter. We need to limit
        % those factors, in case they become unbounded.
        lambda_max = min(supportFunc_(P_iter,x_iter_unit,'upper'), 10);
        lambda_min = -min(supportFunc_(P_iter,-x_iter_unit,'upper'), 10);
        % We now take the midpoint between the two points on the boundary
        % that lie on the line induced by x_iter. That way, we can make
        % sure that the point does not lie on a vertex.
        c_iter = x_iter_unit * (lambda_max+lambda_min)/2;
        
        if any(c_iter)
            % We now re-center P_iter, so that we limit the chances of the
            % origin being on one of the vertices of P_iter;
            P_iter = P_iter - c_iter;
            % We also need to translate all the vectors we have found so
            % far, if we are to find an orthogonal set of vectors:
            if ~isempty(subspace)
                subspace = subspace - c_iter;
            end
        end
        
        % Now, it suffices to add x_iter to the list
        subspace = [subspace,x_iter_unit];
        
    end
    
    % Now that we have constructed the subspace we need, it is time to
    % check what dimension it has; only if it is full dimensional, is the
    % polytope non-degenerate
    k = size(subspace,2);
    if k == dim(P)
        res = true;
        P.fullDim.val = true;
        % If that is the case, we can assume that subspace is the entire
        % space, and so a much easier ONB would be the following
        subspace = eye(dim(P));
    else
        res = false;
        P.fullDim.val = false;
        % It remains to transform this into a ONB
        [Q,R] = qr(subspace);
        subspace = Q(:,1:k);
    end
    
end

% set property
P.fullDim.val = res;

end


% Auxiliary functions -----------------------------------------------------

function [res,subspace] = aux_1D(P)

    % compute vertices
    V = vertices_(P,'lcon2vert');
    % save vertices
    P.V.val = V;
    
    % go over cases...
    if isempty(V)
        % empty set
        res = false; subspace = [];
        P.emptySet.val = true;
        P.fullDim.val = false;
        P.bounded.val = true;
    elseif size(V,2) == 1
        % single vertex
        res = false; subspace = [];
        P.emptySet.val = false;
        P.fullDim.val = false;
        P.bounded.val = true;
    elseif any(isinf(V))
        % unbounded
        res = true; subspace = 1;
        P.emptySet.val = false;
        P.fullDim.val = true;
        P.bounded.val = false;
    else
        % bounded
        res = true; subspace = 1;
        P.emptySet.val = false;
        P.fullDim.val = true;
        P.bounded.val = true;
    end

end


function x = aux_maxNormPerpendicularPolytope(P,X)
    % For a list of vectors X = [x_1,...,x_k], solve the linear program
    %   max ||x||_oo
    %   s.t. x \in P,
    %        forall i: x_i'*x = 0
    %
    % This is equivalent to
    %   max_y max_x y'*x,
    %   s.t. x \in P,
    %        forall i: x_i'*x = 0,
    %        and where y is iterated over all standard vectors +-e_i.
    
    % read out dimension of polytope
    n = dim(P);
    % store maximum objective value and maximizer
    maximum = 0;
    maximizer = [];

    % add constraints forall i: x_i'*x=0
    Aeq = [P.Ae; X'];
    beq = [P.be; zeros(size(X,2),1)];

    P_ = polytope(P.A,P.b,Aeq,beq);

    % loop over all dimensions (for y)
    for i=1:n
        % compute maximum for y = +e_i
        y = zeros(n,1);
        y(i) = 1;
        [res, x] = aux_firstMaximum(P_,y);

        % save maximizer if objective value is larger
        if -res > maximum
            maximum = -res;
            maximizer = x;
        end
        
        % compute maximum for y = -e_i
        y = zeros(n,1);
        y(i) = -1;
        [res, x] = aux_firstMaximum(P_,y);

        % save maximizer if objective value is larger
        if -res > maximum
            maximum = -res;
            maximizer = x;
        end
    end

    % return maximizer
    x = maximizer;
end

function [res, x] = aux_firstMaximum(P,y)

    % linear program options
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    [x,res,exitflag] = linprog(-y, P.A, P.b, P.Ae, P.be, [], [], options);
    
    % If the problem is unbounded, we need to add a constraint, e.g.,
    %   y'*x = 1
    % in order to find a good direction for x.
    if exitflag == -3
        Aeq = [P.Ae; y'];
        beq = [P.be; 1];
        [x,~] = linprog(-y, P.A, P.b, Aeq, beq, [], [], options);
        % set the objective value manually to -Inf to force updating the
        % maximizer
        res = -Inf;
    end
end

% ------------------------------ END OF CODE ------------------------------
