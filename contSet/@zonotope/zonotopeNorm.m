function [res, minimizer] = zonotopeNorm(Z,p)
% zonotopeNorm - computes the norm of the point p w.r.t. the zonotope-norm
%    induced by the zonotope Z (see [1, Definition 4]).
%
% Syntax:
%    res = zonotopeNorm(Z,p)
%
% Inputs:
%    Z - zonotope
%    p - nx1-array, with n the dimension of Z
%
% Outputs:
%    res - zonotope-norm of the point p
%    minimizer - (optional) returns a solution x s.t. Gx = p and for
%                 which norm(x,inf) = zonotopeNorm(Z,p)
%
% Example:
%    c = [0;0];
%    G = [[2 3 0];[2 0 3]];
%    Z = zonotope(c, G);
%    
%    p = rand([2 1]);
%
%    % Set of points that have the same distance to the origin as p, with
%    % respect to the zonotope norm of Z    
%    d = zonotopeNorm(Z, p);
%
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(d*Z,[1,2],'g');
%    plot(p(1),p(2),'rx');
%
% References:
%    [1] A. Kulmburg, M. Althoff. "On the co-NP-Completeness of the
%        Zonotope Containment Problem", European Journal of Control 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       14-May-2021
% Last update:   16-January-2024 (MW, handle edge cases)
%                28-March-2025 (TL, return minimizer)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
narginchk(2,2);

% check input arguments
inputArgsCheck({{Z,'att','zonotope'};
                {p,'att','numeric'}});

% check dimension of Z and p
equalDimCheck(Z, p);

% empty set
if representsa_(Z,'emptySet',eps)
    if isempty(p)
        res = 0; minimizer = []; return
    else
        res = Inf; minimizer = []; return
    end
end

% Retrieve generator-representation of Z
if isempty(Z.G)
    if ~any(p)
        res = 0; minimizer = []; return
    else
        res = Inf; minimizer = []; return
    end
end

% Retrieve dimensions of the generator matrix of Z
[n,numGen] = size(Z.G);

% Set up objective and constraints of the linear program as defined in
% [2, Equation (8)]
problem.f = [1; zeros(numGen,1)];

problem.Aeq = [zeros(n,1), Z.G];
problem.beq = p;

problem.Aineq = [-ones(numGen,1),  speye(numGen); ...
                 -ones(numGen,1), -speye(numGen)];
problem.bineq = zeros(2*numGen,1);

% bounds integrated in inequality constraints
problem.lb = [];
problem.ub = [];

% Solve the linear program: If the problem is infeasible, this means that
% the zonotope must be degenerate, and that the point can not be realized
% as a linear combination of the generators of the zonotope. In that case,
% the norm is defined as Inf. The same goes when the problem is 'unbounded'
[minimizer_p,res,exitflag] = CORAlinprog(problem);
if exitflag == -2 || exitflag == -3
    res = Inf; minimizer = []; return;
elseif exitflag ~= 1
    % In case anything else went wrong, throw out an error
    throw(CORAerror('CORA:solverIssue'));
end

% all good, compute minimizer
    
% remove first entry from minimizer as obtained by lp
minimizer = minimizer_p(2:end);


% ------------------------------ END OF CODE ------------------------------
