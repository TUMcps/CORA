function [res,S] = representsa_(cZ,type,tol,varargin)
% representsa_ - checks if a constrained zonotope can also be represented
%    by a different set, e.g., a special case
%
% Syntax:
%    res = representsa_(cZ,type,tol)
%    [res,S] = representsa_(cZ,type,tol)
%
% Inputs:
%    cZ - conZonotope object
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

% Authors:       Mark Wetzlinger, Niklas Kochdumper
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(cZ,type);
else
    [empty,res,S] = representsa_emptyObject(cZ,type);
end
if empty; return; end

% dimension
n = dim(cZ);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        res = representsa_(interval(cZ),'origin',tol);
        if nargout == 2 && res
            S = zeros(n,1);
        end

    case 'point'
        res = representsa_(interval(cZ),'point',tol);
        if nargout == 2 && res
            S = center(cZ);
        end

    case 'capsule'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));

    case 'conHyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));
        
    case 'conPolyZono'
        % obviously true
        res = true;
        if nargout == 2
            S = conPolyZono(cZ);
        end

    case 'conZonotope'
        % obviously true
        res = true;
        if nargout == 2
            S = cZ;
        end

    case 'ellipsoid'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));

    case 'halfspace'
        % constrained zonotopes cannot be unbounded
        res = false;

    case 'interval'
        res = isempty(cZ.A) && representsa_(zonotope(cZ),'interval',tol);
        if nargout == 2 && res
            S = interval(cZ);
        end

    case 'levelSet'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));

    case 'polytope'
        res = true;
%         if nargout == 2 && res
%             S = polytope(cZ);
%         end

    case 'polyZonotope'
        res = true;
        if nargout == 2
            S = polyZonotope(cZ);
        end

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        res = true;
        if nargout == 2
            S = zonoBundle(cZ);
        end

    case 'zonotope'
        res = isempty(cZ.A);
        % note: there may be cases with constraints that can still be
        % represented by zonotopes
        if nargout == 2 && res
            S = zonotope(cZ);
        end

    case 'hyperplane'
        % constrained zonotopes cannot be unbounded (unless 1D, where
        % hyperplane is also bounded)
        res = n == 1;

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conZonotope to ' type ' not supported.']));

    case 'convexSet'
        res = true;

    case 'emptySet'
        res = aux_isEmptySet(cZ,tol);
        if nargout == 2 && res
            S = emptySet(n);
        end

    case 'fullspace'
        % constrained zonotopes cannot be unbounded
        res = false;

end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isEmptySet(cZ,tol)

% check if the (in general, enclosing) zonotope is already empty
if representsa_(zonotope(cZ.c,cZ.G),'emptySet',tol)
    res = true; return
end

% if there are no constraints, we are finished
if isempty(cZ.A)
    res = false; return
end

% approach: if the constraints are satisfiable, there is at least one value
% for the beta factors and, thus, the constrained zonotope is non-empty
% (note: one should be able to replace this by checking
%   P = interval(-ones(nrGen,1),ones(nrGen,1)) & polytope([],[],cZ.A,cZ.b);
%   res = representsa(P,'emptySet',tol);
% ...do this diligently sometime in the future)

% functions below do not support sparse matrices
if issparse(cZ.A)
    cZ.A = full(cZ.A);
end

% null space of the constraints
Neq=null(cZ.A);

% find a single point that satisfies the constraints
x0 = pinv(cZ.A)*cZ.b;

if norm(cZ.A*x0 - cZ.b) > 1e-10*norm(cZ.b)  % infeasible
    % note: the tolerance above must be hardcoded to some non-zero value
    res = true;
    return
end

if isempty(Neq)  % null space empty -> constraints admit a single point
    % check if the single point for beta satisfies the unit cube
    nrGen = size(cZ.G,2);
    res = ~contains_(interval(-ones(nrGen,1),ones(nrGen,1)),x0,'exact',tol,0,false,false);
    return
end

% check if the null-space intersects the unit-cube
[nrCon,nrGen] = size(cZ.A);
unitCube = interval(-ones(nrGen,1),ones(nrGen,1));

% loop over all constraints (= hyperplanes)
for i = 1:nrCon
    % hyperplane from a constraint does not intersect the unit cube
    % -> set is empty
    if ~isIntersecting_(polytope(cZ.A(i,:),cZ.b(i)),unitCube,'exact',tol)
        res = true; return
    end
end

% use linear programming to check if the constrained zonotope is
% empty (this seems to be more robust than the previous solution
% using the polytope/isempty function)
if nrCon >= 1

    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off', ...
                                'OptimalityTolerance',1e-10);
    end

    problem.f = ones(nrGen,1);
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = cZ.A;
    problem.beq = cZ.b;
    problem.ub = ones(nrGen,1);
    problem.lb = -problem.ub;
    problem.options = options;
    
    [~,~,exitflag] = CORAlinprog(problem);
    
    res = exitflag == -2;
end

end

% ------------------------------ END OF CODE ------------------------------
