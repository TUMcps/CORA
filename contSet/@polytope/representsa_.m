function [res,S] = representsa_(P,type,tol,varargin)
% representsa_ - checks if a polytope can also be represented by a
%    different set, e.g., a special case
%
% Syntax:
%    res = representsa_(P,type,tol)
%    [res,S] = representsa_(P,type,tol)
%
% Inputs:
%    P - polytope object
%    type - other set representation or 'origin', 'point', 'hyperplane'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    S - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa

% Authors:       Mark Wetzlinger, Niklas Kochdumper, Victor Gassmann
% Written:       25-July-2023
% Last update:   01-August-2023 (MW, support fullspace comparison)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [emptyObj,res] = representsa_emptyObject(P,type);
else
    [emptyObj,res,S] = representsa_emptyObject(P,type);
end
if emptyObj; return; end

% dimension
n = dim(P);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        % quick check: is origin contained?
        res = false;
        if contains_(P,zeros(n,1),'exact',tol)
            % definitely not empty if the origin is contained
            P.emptySet.val = false;
            % check if only origin contained
            res = norm_(interval(P),2) <= tol;
            % set is degenerate if it's only the origin
            if res
                P.fullDim.val = false;
            end
            if nargout == 2 && res
                S = zeros(n,1);
            end
        end        

    case 'point'
        if n == 1
            V = vertices(P);
            res = size(V,2) == 1;
        else
            [fulldim,subspace] = isFullDim(P);
            res = ~fulldim && isempty(subspace);
        end
        % set is degenerate if it's only a single point
        if res
            P.fullDim.val = false;
        end
        if nargout == 2 && res
            % only one point
            S = center(P);
        end

    case 'capsule'
        % true if 1D and bounded
        % note: also true if polytope is a bounded line
        res = n == 1 && isBounded(P);
        if nargout == 2
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from polytope to ' type ' not supported.']));
        end

    case 'conHyperplane'
        [res,S] = aux_isConHyperplane(P);

    case 'conPolyZono'
        res = isBounded(P);
        if nargout == 2 && res
            S = conPolyZono(P);
        end

    case 'conZonotope'
        % always true for bounded polytopes
        res = isBounded(P);
        if nargout == 2 && res
            S = conZonotope(P);
        end

    case 'ellipsoid'
        % only an ellipsoid if 1D and bounded or a single point
        res = (n == 1 && isBounded(P) && ~representsa_(P,'emptySet',tol)) ...
            || representsa_(P,'point',tol);
        if nargout == 2 && res
             S = ellipsoid(P);
        end

    case 'halfspace'
        % if only a single irredundant inequality constraint given
        P = compact_(P,'all',tol);
        res = size(P.A,1) == 1 && isempty(P.Ae);
        if nargout == 2 && res
            S = halfspace(P.A,P.b);
        end

    case 'interval'
        if nargout == 1
            res = aux_isInterval(P,tol);
        elseif nargout == 2
            [res,S] = aux_isInterval(P,tol);
        end

    case 'levelSet'
        res = true;
        if nargout == 2
            S = levelSet(P);
        end

    case 'polytope'
        % obviously true
        res = true;
        if nargout == 2
            S = P;
        end

    case 'polyZonotope'
        % only true if polytope is bounded
        res = isBounded(P);
        if nargout == 2 && res
            S = polyZonotope(P);
        end

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        % zonotope bundles can represent any bounded convex set
        res = isBounded(P);
        if nargout == 2 && res
            S = zonoBundle(P);
        end

    case 'zonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of polytope to ' type ' not supported.']));

    case 'hyperplane'
        res = isempty(P.A) && size(P.Ae,1) == 1;
        if res
            % hyperplanes are unbounded and non-empty
            P.bounded.val = false;
            P.emptySet.val = false;
        end
        if nargout == 2 && res
            S = conHyperplane(P.Ae,P.be);
        end

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of polytope to ' type ' not supported.']));

    case 'emptySet'
        res = aux_isEmptySet(P);
        P.emptySet.val = res;
        if res
            P.bounded.val = true;
            P.fullDim.val = false;
            P.V.val = zeros(n,0);
            P.minVRep.val = true;
        end
        if nargout == 2 && res
            S = emptySet(n);
        end

    case 'fullspace'
        % all constraints must be trivially fulfilled, that is,
        %    A * x <= a, a >= 0, and  Ae* x = 0
        res = all(all(withinTol(P.A,0,tol))) ...
            && all(P.b > 0 | withinTol(P.b,0,tol)) ...
            && all(all(withinTol(P.Ae,0,tol))) ...
            && all(withinTol(P.be,0,tol));
        % fullspaces are always unbounded, non-empty and full-dimensional
        if res
            P.bounded.val = false;
            P.emptySet.val = false;
            P.fullDim.val = true;
        end
        if nargout == 2 && res
            S = fullspace(n);
        end

end

end


% Auxiliary functions -----------------------------------------------------

function [res,p] = aux_isZero(P,varargin)
% former isZero function...

    [tol,mode] = setDefaultValues({0,'approx'},varargin);
    
    % parse input arguments
    inputArgsCheck({{P,'att','polytope'};
                    {tol,'att','numeric',{'nonnegative','scalar','nonnan'}}; ...
                    {mode,'str',{'approx','approx:para','exact'}}});
    
    % default result (if res = false)
    p = [];

    % empty set case
    if isempty(P)
        res = false; 
        return;
    end
    
    if strcmp(mode,'approx')
        % requires 2*n LPs
        r = norm(interval(P));
    
    elseif strcmp(mode,'approx:para')
        % compute parallelotope enclosure of P (requires SDP)
        Zp = zonotope(P);
        
        % radius of parallelotope is 
        r = norm(center(Zp)+sum(abs(generators(Zp)),2));
    
    elseif strcmp(mode,'exact')
        throw(CORAerror("CORA:noExactAlg",P));
    end
    
    res = r <= tol;
    if res
        p = zeros(n,1);
    end

end

function res = aux_isEmptySet(P)
% checks if a polytope is empty
%       A*x <= b
%    is infeasible if
%        min  b'*y
%        s.t. A'*y =  0,
%                y >= 0
%    is feasible. If problem is either unbounded (below since minimization)
%    or b'*y < 0, the polytope is empty. If the problem is infeasible or
%    b'*y == 0, the polytope is not empty.

    % check if emptiness already known
    if ~isempty(P.emptySet.val)
        res = P.emptySet.val; return
    end
    
    % if polytope has V-Rep, then not empty
    if ~isempty(P.V.val)
        res = false; return
    end
    
    % dimension
    n = dim(P);
    
    % special case: 1D
    if n == 1    
        % remove redundancies (fast for 1D), this also finds out whether
        % the set is empty or not
        P_min = compact_(P,'all',1e-9);
        res = P_min.emptySet.val;
        return
    end
    
    % quick check: are there constraints of the form 0*x <= b with b < 0?
    zero_rows = find(all(withinTol(P.A,0),2));
    if any(zero_rows) && any(P.b(zero_rows) < 0)
        res = true; return
    end
    
    
    % number of inequality and equality constraints
    nrConIneq = size(P.A,1);
    nrConEq = size(P.Ae,1);
    
    % normalize constraints
    P = normalizeConstraints(P,'A');
    
    % read out properties
    A = P.A; b = P.b;
    Ae = P.Ae; be = P.be;
    
    if ~isempty(Ae)
        % polytope can be already be declared empty if there is not point
        % that satisfies the set of equality constraints
        if rank(Ae) == min([nrConEq,n])
            % if the rank of the equality constraint matrix is smaller than
            % the minimum of the number of equality constraints and the
            % dimension, then there are still degrees of freedom
            
            % check for aligned vectors: pre-compute dot product
            dotprod_norm = Ae * Ae';
            
            for i=1:nrConEq
                % if there are aligned constraints, then they have to have
                % the same value in be (all but one of them are redundant),
                % otherwise the polytope is the empty set because there is
                % no point that can satisfy to parallel Ae*x = be at the
                % same time
    
                % check for aligned vectors
                alignedConstraints = all(...
                    withinTol(Ae(i,:) - dotprod_norm(:,i)*Ae(i,:),...
                    zeros(1,n)), 2);
            
                if nnz(alignedConstraints) > 1
                    % at least two constraints are aligned
                    if ~all(withinTol(be(alignedConstraints),be(i)))
                        % polytope is the empty set
                        res = true;
                        P.emptySet.val = true;
                        return
                    end
                end
            end
    
        end
    
        % rewrite equality constraints as inequality constraints
        A = [A; Ae; -Ae];
        b = [b; be; -be];
        % update number of inequality constraints
        nrConIneq = length(b);
    end
    
    
    % solve the dual problem using linear programming
    
    if isSolverInstalled('mosek')
    
        % call MOSEK
        res = msklpopt(b,[-eye(nrConIneq);A';-A'],...
            -Inf(nrConIneq+2*n,1),zeros(nrConIneq+2*n,1),...
            [],[],[],'minimize echo(0)');
    
        % read out solution
        if strcmp(res.sol.itr.prosta,'DUAL_INFEASIBLE')
            % unbounded (from below since minimization)  =>  empty
            res = true;
        elseif strcmp(res.sol.itr.prosta,'PRIMAL_INFEASIBLE')
            % no point could be found -> not empty
            res = true;
        elseif strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
            if res.sol.itr.dobjval > 0 || withinTol(res.sol.itr.dobjval,0)
                % optimal objective value is >= 0  =>  not empty
                res = false;
            else
                % optimal objective value is < 0  =>  empty
                res = true;
            end
        end
    else
    
        % set options for MATLAB linprog
        persistent options
        if isempty(options)
            options = optimoptions('linprog','display','off', ...
                               'OptimalityTolerance',1e-10);
        end
        
        % solve linear program (all inequalities)
        [x,~,exitflag] = linprog(b',...
            [-eye(nrConIneq);A';-A'],zeros(nrConIneq+2*n,1),...
            [],[],[],[],options);
        
        % read out solution
        if exitflag == -2 || (exitflag > 0 && b'*x >= -1e-10)
            % if this problem is infeasible or if optimal objective value
            % is 0, the polytope is not empty
            res = false;
        elseif exitflag == -3 || (exitflag > 0 && b'*x < 1e-10)
            % if problem is unbounded (below since minimization) or
            % objective value is smaller zero, polytope is empty
            res = true;
        else
            throw(CORAerror('CORA:solverIssue','linprog'))
        end
    end

end

function [res,I] = aux_isInterval(P,tol)

    % assume false
    res = false; I = [];

    % TODO: rewrite this function
    
    % fast initial check
    if ~all(sum(abs(sign(P.A)),2) == 1)
        return;
    end
    
    % construct equivalent interval
    n = dim(P); lb = -Inf(n,1); ub = Inf(n,1);
    
    % loop over each dimension
    for i = 1:n
        % find halfspaces along this dimension
        ind = find(P.A(:,i) ~= 0);
        % loop over each halfspace
        for j = 1:length(ind)
            if P.A(ind(j),i) > 0
                % upper bound in this dimension
                ub(i) = min(ub(i),P.b(ind(j))/P.A(ind(j),i)); 
            else
                % lower bound in this dimension
                lb(i) = max(lb(i),P.b(ind(j))/P.A(ind(j),i));
            end
        end
    end
    
    % check if the interval is bounded
    if nargout == 2
        I = interval(lb,ub); res = true; 
    end

end

function [res,hyp] = aux_isConHyperplane(P)

    % assume false
    res = false;
    hyp = [];

    % check if equality constraints exist
    if ~isempty(P.Ae)
        
        % first equality constraint defines hyperplane
        c = P.Ae(1,:);
        d = P.be(1);
        
        % convert all other equality constraints to inequality constraints
        A = [P.A; P.Ae(2:end,:); -P.Ae(2:end,:)];
        b = [P.b; P.be(2:end); -P.be(2:end)];
        
        % construct constrained hyperplane
        res = true;
        if isempty(A) && isempty(b)
            hyp = conHyperplane(c,d);
        else
            hyp = conHyperplane(c,d,A,b);
        end
        
        
    else
        
        % remove redundant halfspaces
        P = compact_(P,'all',1e-9);
        
        % normalize the inequality constraints
        A = P.A;
        b = P.b;
        
        n = sum(A.^2,2);
        A = diag(1./n)*A;
        b = b./n;
        
        % check for the case A*x <= b && A*x >= b => A*x == b
        for i = 1:size(A,1)-1
            
            % check if normal vectors are identical
            [result,idx] = ismembertol(-A(i,:),A(i+1:end,:),'ByRows',true);
            idx = idx(1) + i;
         
            if result
                
                % check if offsets are identical
                if b(i) == -b(idx)
                   
                    % equality constraint detected -> construct hyperplane
                    c = A(i,:);
                    d = b(i);
                    
                    A([i,idx],:) = [];
                    b([i,idx],:) = [];
                    
                    res = true;
                    if isempty(A) && isempty(b)
                        hyp = conHyperplane(c,d);
                    else
                        hyp = conHyperplane(c,d,A,b);
                    end
                    return;
                end
            end
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
