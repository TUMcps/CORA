function res = testLong_polytope_polytope
% testLong_polytope_polytope - unit test function of constructor
%
% Syntax:
%    res = testLong_polytope_polytope
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       05-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 50;

for i=1:nrTests

    % random dimension
    n = randi(10);

    % sample semi-random vertices (no redundancies)
    V = -1+2*rand(n,2*n);
    V = V + [diag(10*ones(n,1)), -diag(10*ones(n,1))];

    % compute H-representation
    P = polytope(V);

    % vertices must be contained in polytope
    if ~all(contains(P,V,'exact',1e-14))
        res = false; return
    end

    % sample vertices (likely with redundancies)
    A = rand(n,n);
    c = randn(n,1);
    sigma = 0.5*(A+A') + n*eye(n);
    V = mvnrnd(c,sigma,2*n)';

    % compute H-representation
    P = polytope(V);

    % vertices must be contained in polytope
    if ~all(contains(P,V,'exact',1e-14))
        res = false; return
    end

end

% ------------------------------ END OF CODE ------------------------------
