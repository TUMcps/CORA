function c = center(P)
% center - Computes the Chebyshev center of a given polytope
%
% Syntax:
%    c = center(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    c - Chebyshev center of the polytope
%
% Example:
%    P = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    c = center(P);
%    
% Reference: MPT-Toolbox https://www.mpt3.org/    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       28-March-2022
% Last update:   14-December-2022 (MW, add call to MOSEK)
%                27-July-2023 (MW, add fast method for 1D)
%                02-January-2024 (MW, fix fully empty polytopes)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimension
n = dim(P);

% Check if polytope is empty
if representsa_(P,'fullspace',0)
    c = NaN(n,1); return;
end

% fast and simple computation for 1D
if n == 1
    % compute vertices
    V = vertices_(P,'lcon2vert');
    if isempty(V)
        % empty set
        c = [];
    elseif size(V,2) == 1
        % only one vertex -> center
        c = V;
    elseif any(isinf(V))
        % unbounded
        c = NaN;
    else
        % bounded, two vertices
        c = mean(V);
    end
    return
end

% check whether there are only equalities
if isempty(P.A) && ~isempty(P.Ae)
    % three outcomes: unbounded, single point, infeasible

    % minimal halfspace representation: if two constraints are aligned and
    % cannot be fulfilled at the same time, an empty polytope is returned
    P_ = compact_(P,'Ae',1e-10);
    if ~isempty(P.emptySet.val) && P.emptySet.val
        c = double.empty(n,0); return;
    end

    % all constraints now are linearly independent, hence the relation of 
    % system dimension and number of constraints determines the solution
    if size(P_.Ae,1) < n
        % underdetermined -> unbounded
        c = NaN(n,1);
    elseif size(P_.Ae,1) > n
        % overdetermined -> no solution
        c = double.empty(n,0);
    else
        % same number of constraints as system dimension -> single point
        c = P_.Ae \ P_.be;
    end
    return
end

% dimension and number of (in)equalities
nrEq = size(P.Ae,1);

% 2-Norm of each row
A_norm = sqrt(sum(P.A.^2,2));

% extend inequality and equality constraints by one column
A = [P.A A_norm];
Ae = [P.Ae zeros(nrEq,1)];

% cost function for linear program: minimize 2-norm of constraints
f = [zeros(n,1); -1];

% different solvers
if isSolverInstalled('mosek')
    
    % rewrite for MOSEK syntax
    a = [A; Ae];
    blc = [-Inf(size(A,1),1); P.be];
    buc = [P.b; P.be];
    blx = [-Inf(n,1); 0];
    bux = [];

    % call MOSEK
    
    res = msklpopt(f,a,blc,buc,blx,bux,[],'minimize echo(0)');

    % read out solution
    if strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
        c = res.sol.itr.xx(1:n);
    elseif strcmp(res.sol.itr.prosta,'PRIMAL_INFEASIBLE')
        % set is empty
        c = [];
    elseif strcmp(res.sol.itr.prosta,'DUAL_INFEASIBLE')
        % unbounded
        c = NaN(n,1);
    end

else
    % MATLAB linprog

    % linear program options
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    
    % Solve Linear Program
    [c,val,exitflag] = linprog(f,A,P.b,Ae,P.be,[-Inf(n,1);0],[],options);

    if exitflag == 1
        % truncate solution
        c = c(1:n);
    elseif exitflag == -2
        % set is empty
        c = double.empty(n,0);
    elseif exitflag == -3
        % unbounded
        c = NaN(n,1);
    elseif exitflag < 0
        throw(CORAerror('CORA:solverIssue'));
    end

end

% ------------------------------ END OF CODE ------------------------------
