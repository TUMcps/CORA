function [res, subspace] = isFullDim(P)
% isFullDim - check if a polytope is full-dimensional
%
% Syntax:  
%    res = isFullDim(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    res - true/false
%    subspace - (optional) Returns a set of orthogonal unit vectors
%               x_1,...,x_k such that P is strictly contained in
%               center(P)+span(x_1,...,x_k) (here, 'strictly' means that k
%               is minimal).
%               Note that if P is just a point, subspace=[].
%
% Example:
%    P1 = mptPolytope([1 0;0 1;-1 0; 0 -1],[1;1;0;0]);
%    P2 = mptPolytope([1 0;0 1;-1 0; 0 -1],[1;0;0;0]);
%
%    isFullDim(P1)
%    isFullDim(P2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Author:       Niklas Kochdumper, Viktor Kotsev, Mark Wetzlinger, Adrian Kulmburg
% Written:      02-January-2020 
% Last update:  08-August-2022
%               05-December-2022 (VK, new method)
%               14-December-2022 (MW, add call to MOSEK)
%               25-May-2023 (AK, added method to compute subspace)
%               28-May-2023 (MW, add quick check for pairwise annihilation)
% Last revision:---

%------------- BEGIN CODE --------------

% too many input arguments
if nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% check if the polytope is completely empty
if isemptyobject(P)
    res = false;
    subspace = [];
    return
end

if nargout < 2
    % If the user only wants to know whether the polytope is degenerate, we
    % can check this rapidly using linear programming

    % Existence of equality constraints implies degenerate case
    if ~isempty(P.P.Ae)
        res = false;
        return
    end

    % number of inequality constraints and dimension
    [nrCon,n] = size(P.P.A);

    % normalize constraints
    P = normalizeConstraints(P,'A');

    % quick check: any pair of inequalities of type
    %   ax <= b, ax >= b_ (-ax <= -b_), with b_ >= b?
    dotprod = tril(P.P.A * P.P.A');
    idxAntiParallel = find(withinTol(dotprod,-1));
    for i=1:length(idxAntiParallel)
        % read out matrix indexing from linear indexing
        [r,c] = ind2sub([nrCon,nrCon],idxAntiParallel(i));
        % check offset values, two perspectives:
        %  ax <= b1, -ax <= b2 <=> ax >= -b2  ->  -b2 >= b1?
        % -ax <= b1 <=> ax >= -b1, ax <= b2   ->  -b1 >= b2 <=> -b2 >= b1?
        % ...the one with the reversed sign has to be larger (or equal)!
        if -P.P.b(c) > P.P.b(r) || withinTol(-P.P.b(c),P.P.b(r))
            res = false;
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
    A = [P.P.A A_norm];

    % cost function for linear program: minimize 2-norm of constraints
    f = [zeros(n,1); -1];

    % different solvers
    if isMosek

        % call MOSEK
        resMOSEK = msklpopt(f,A,-Inf(nrCon,1),P.P.b,...
            [],[],[],'minimize echo(0)');

        % read out solution
        if strcmp(resMOSEK.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
            % check if point is strictly interior (radius of center = 0)
            if withinTol(0,resMOSEK.sol.itr.dobjval,1e-9)
                res = false;
            else
                res = true;
            end
        elseif strcmp(resMOSEK.sol.itr.prosta,'PRIMAL_INFEASIBLE')
            % set is empty
            res = false;
        elseif strcmp(resMOSEK.sol.itr.prosta,'DUAL_INFEASIBLE')
            % numerical issues seem to come up only in the case of
            % full-dimensional unbounded polytopes
            res = true;
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
        [~,r,exitflag] = linprog(f,A,P.P.b,[],[],[-Inf(n,1);0],[],options);

        if exitflag == 1
            % check if point is strictly interior (radius of center = 0)
            if isempty(r)
                res = false;
            elseif withinTol(0,r,1e-9)
                res = false;
            else
                res = true;
            end
        elseif exitflag == -2
            % set is empty
            res = false;
        elseif exitflag == -3
            % numerical issues seem to come up only in the case of
            % full-dimensional unbounded polytopes
            res = true;
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
        res = 0;
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
        % If the solution is x_iter=0, then this means that there are no
        % non-trivial vectors that are perpendicular to those we have
        % found, and so we can stop the iteration
        if norm(x_iter,Inf) <= 1e-5
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
        % If that is the case, we can assume that subspace is the entire
        % space, and so a much easier ONB would be the following
        subspace = eye(dim(P));
    else
        res = false;
        % It remains to transform this into a ONB
        [Q,R] = qr(subspace);
        subspace = Q(:,1:k);
    end
    
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
    Aeq = [P.P.Ae; X'];
    beq = [P.P.be; zeros(size(X,2),1)];

    polyStruct = struct('A',P.P.A,'b',P.P.b,'Ae',Aeq,'be',beq);
    P_ = mptPolytope(polyStruct);

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
    [x,res,exitflag] = linprog(-y, P.P.A, P.P.b, P.P.Ae, P.P.be, [], [], options);
    
    % If the problem is unbounded, we need to add a constraint, e.g.,
    %   y'*x = 1
    % in order to find a good direction for x.
    if exitflag == -3
        Aeq = [P.P.Ae; y'];
        beq = [P.P.be; 1];
        [x,~] = linprog(-y, P.P.A, P.P.b, Aeq, beq, [], [], options);
        % set the objective value manually to -Inf to force updating the
        % maximizer
        res = -Inf;
    end
end

%------------- END OF CODE --------------