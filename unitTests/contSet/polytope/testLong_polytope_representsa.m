function res = testLong_polytope_representsa
% testLong_polytope_representsa - unit test function for empty check
%
% Syntax:
%    res = testLong_polytope_representsa()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       29-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 100;

for i=1:nrTests
    
    % random dimension
    n = randi(10);

    % create interval-like polytopes: halfspaces almost axis-aligned
    % (fill zero values with very small factors)

    % random matrices
    M1 = 0.001*rand(n);
    M2 = 0.001*rand(n);
    A = [M1 - diag(diag(M1)) + eye(n); M2 - diag(diag(M2)) - eye(n)];

    % offset vector large enough so that result will not be empty
    b = 1000*ones(2*n,1);

    % init polytope
    P = polytope(A,b);

    % emptiness check
    if representsa(P,'emptySet')
        res = false; return
    end

    % init polytope where one pair of constraints is unfulfillable

    % number of constraints
    nrCon = randi([2,2*n]);

    % random dimension where constaint is unfulfillable
    randDim = randi([1,n]);

    % init pair of constraints
    I = eye(n);
    A = [I(:,randDim), -I(:,randDim)];
    b = [-rand; -rand];

    % fill up remaining constraints
    A = [A, zeros(n,nrCon-2)];
    b = [b; zeros(nrCon-2,1)];
    for j=1:nrCon-2
        % all entries in A between 0 and 0.01
        A(:,j+2) = 0.01*rand(n,1);
        % all b values between 500 and 1000
        b(j+2) = 500 + 500*rand(1);
    end

    % init polytope
    P = polytope(A',b);

    % emptiness check
    if ~representsa(P,'emptySet')
        res = false; return
    end

end

% ------------------------------ END OF CODE ------------------------------
