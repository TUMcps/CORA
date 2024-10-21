function D_ksi = priv_ksi_optimizer(cZ)    
% priv_ksi_optimizer - determine the tightened domains for the zonotope
%    factors ksi by solving a linear program
%
% Syntax:
%    Dksi = priv_ksi_optimizer(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    D_ksi - new tightened domains for the zonotope factors (class: interval)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:       11-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(cZ.A)
    % no constraints -> return unit-cube as domain for ksi
    n = size(cZ.G,2);
    D_ksi = interval(-ones(n,1),ones(n,1));
    return
end

% object properties  
n = size(cZ.A, 2);
A = cZ.A;
b = cZ.b;

% ksi in [-1, 1]
lb = -ones(n,1);
ub = ones(n,1);

% initialize ksi
ksi_min = zeros(n,n);
ksi_max = zeros(n,n);

% init linprog struct
problem.Aineq = [];
problem.bineq = [];
problem.Aeq = A;
problem.beq = b;
problem.lb = lb;
problem.ub = ub;

for i = 1:n

    % min/max ksi_i
    f = zeros(n, 1);
    f(i, 1) = 1;

    % minimize (Equation (25) in [1])
    problem.f = f;
    [x, ~, flag_min] = CORAlinprog(problem);
    ksi_min(:,i) = x;

    % maximize  (Equation (26) in [1])
    problem.f = -f;
    [x, ~, flag_max] = CORAlinprog(problem);
    ksi_max(:,i) = x;
    if flag_min ~= 1 || flag_max ~= 1
        throw(CORAerror('CORA:solverIssue'));
    end
end

% delete duplicates
ksi = unique([ksi_min, ksi_max]','rows')';

% calculate tightened domain for the zonotope factors
D_ksi = interval(min(ksi,[],2),max(ksi,[],2));

% ------------------------------ END OF CODE ------------------------------
