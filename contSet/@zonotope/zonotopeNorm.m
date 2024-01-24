function res = zonotopeNorm(Z,p)
% zonotopeNorm - computes the norm of the point p w.r.t. the zonotope-norm
%    induced by the zonotope Z (see [1, Definition 4])
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{Z,'att','zonotope'};
                {p,'att','numeric'}});

% empty set
if representsa_(Z,'emptySet',eps)
    if isempty(p)
        res = 0; return
    else
        res = Inf; return
    end
end

% Retrieve generator-representation of Z
G = Z.generators;
if isempty(G)
    if ~any(p)
        res = 0; return
    else
        res = Inf; return
    end
end

% Retrieve dimensions of the generator matrix of Z
n = size(G, 1);
m = size(G, 2);

% Set up objective and constraints of the linear program as defined in
% [2, Equation (8)]
f = [1;zeros([m 1])];

Aeq = [zeros([n 1]) G];
beq = p;

Aineq1 = [-ones([m 1]) eye(m)];
Aineq2 = [-ones([m 1]) -eye(m)];

Aineq = [Aineq1; Aineq2];
bineq = zeros([2*m 1]);

% Suppress solver output
persistent options
if isempty(options)
    options = optimoptions('linprog', 'Display', 'none');
end

% Solve the linear program
[minimizer_p, res, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, [], [], options);

% If the problem is not feasible, this means that the zonotope must be
% degenerate, and that the point can not be realized as a linear
% combination of the generators of the zonotope. In that case, the norm is
% defined as inf
if exitflag == -2
    res = inf;
% In case anything else went wrong, throw out an error
elseif exitflag ~= 1
    throw(CORAerror('CORA:solverIssue'));
end

% ------------------------------ END OF CODE ------------------------------
