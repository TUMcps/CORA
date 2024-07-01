function [res,subspace] = isFullDim(P)
% isFullDim - checks if the dimension of the affine hull of a polytope is
%    equal to the dimension of its ambient space; additionally, one can
%    obtain a basis of the subspace in which the polytope is contained
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
%                03-May-2024 (MW, fix subspace computation for V-rep)
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
    rankV = rank(P.V.val - mean(P.V.val,2),1e-10);
    res = rankV == n;
    if res
        subspace = eye(n);
    else
        [Q,R] = qr(P.V.val - mean(P.V.val,2));
        subspace = Q(:,1:rankV);
    end
    P.fullDim.val = res;
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

        % init linprog struct
        problem.f = f;
        problem.Aineq = A;
        problem.bineq = P.b;
        problem.lb = [-Inf(n,1);0];

        % linear program options
        persistent options
        if isempty(options)
            options = optimoptions('linprog','display','off');
        end

        problem.solver = 'linprog';
        problem.options = options;

        % Solve Linear Program
        [~,r,exitflag] = linprog(problem);

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

    % since we have to evaluate many linear programs anyway, we can quickly
    % check if the polytope is completely empty -> also degenerate
    if representsa_(P,'emptySet',1e-10)
        res = false;
        subspace = [];
        return
    end
    
    % Let's first check whether the polytope is empty, by searching for one
    % point that is contained in P:
    x0 = aux_maxNormPerpendicularPolytope(P,zeros([dim(P) 0]));
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
    subspace = zeros([dim(P) 0]); % Need to have this form to avoid bugs
    
    for i = 1:n
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
        
        x_iter_unit = x_iter / norm(x_iter);
        
        % Now, it suffices to add x_iter to the list
        subspace = [subspace x_iter_unit];
        
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
    %   s.t. x + Xw \in P,
    %        forall i: x_i'*x = 0
    %
    % This is equivalent to
    %   max_y max_x y'*x,
    %   s.t. x + Xw \in P,
    %        forall i: x_i'*x = 0,
    %        and where y is iterated over all standard vectors +-e_i.
    
    % read out dimension of polytope
    n = dim(P);
    % store maximum objective value and maximizer
    maximum = 0;
    maximizer = [];

    % loop over all dimensions (for y)
    for i=1:n
        % compute maximum for y = +e_i
        y = zeros(n,1);
        y(i) = 1;
        [res, x] = aux_firstMaximum(P,y,X);

        % save maximizer if objective value is larger
        if -res > maximum
            maximum = -res;
            maximizer = x;
        end
        
        % compute maximum for y = -e_i
        y = zeros(n,1);
        y(i) = -1;
        [res, x] = aux_firstMaximum(P,y,X);

        % save maximizer if objective value is larger
        if -res > maximum
            maximum = -res;
            maximizer = x;
        end
    end

    % return maximizer
    x = maximizer;
end

function [res, x] = aux_firstMaximum(P,y,X)

    % linear program options
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    
    % Extend to have x + Xw \in P
    P_A_extended = [P.A P.A*X];
    P_Aeq_extended = [P.Ae P.Ae*X];
    y_extended = [y; zeros([size(X,2) 1])];
    
    % add constraints forall i: x_i'*x=0
    P_Aeq_extended = [P_Aeq_extended; X' zeros([size(X,2) size(X,2)])];
    P_beq_extended = [P.be; zeros(size(X,2),1)];

    % init linprog struct
    problem.f = -y_extended;
    problem.Aineq = P_A_extended;
    problem.bineq = P.b;
    problem.Aeq = P_Aeq_extended;
    problem.beq = P_beq_extended;
    problem.solver = 'linprog';
    problem.options = options;
    
    [x,res,exitflag] = linprog(problem);
    
    % If the problem is unbounded, we need to add a constraint, e.g.,
    %   y'*x = 1
    % in order to find a good direction for x.
    if exitflag == -3
        problem.Aeq = [P_Aeq_extended; y_extended'];
        problem.beq = [P_beq_extended; 1];
        [x,~] = linprog(problem);
        % set the objective value manually to -Inf to force updating the
        % maximizer
        res = -Inf;
    end

    x = x(1:dim(P));
end

% ------------------------------ END OF CODE ------------------------------
