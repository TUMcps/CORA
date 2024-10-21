function res = testLong_polytope_isFullDim
% testLong_polytope_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = testLong_polytope_isFullDim
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
% Written:       28-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% number of tests
nrTests = 20;

for i=1:nrTests

    % random dimension
    n = randi(8);

    % init polytope from parallelotope (quick conversion)
    Z = zonotope.generateRandom('Dimension',n,'NrGenerators',n);
    P = polytope(Z);

    % check if zonotope was full-dimensional
    tol = 1e-4;
    if ~isFullDim(Z,tol)
        continue
    end

    % assert non-degeneracy (tol in Z is less strict?)
    assertLoop(isFullDim(P,tol/10),i);

    % add randomly-oriented inequality of type ax <= b, -ax <= -b
    con_ = randn(1,n);
    offset_ = randn(1);
    P_ = polytope([P.A; con_; -con_],[P.b; offset_; -offset_]);

    % assert degeneracy
    assertLoop(~isFullDim(P_),i);

    % lift dimension via equality constraint x_n+1 = 0
    P_ = polytope([P.A,zeros(length(P.b),1)],P.b,[zeros(1,n),1],0);
    [res_,subspace] = isFullDim(P_);
    % has to be degenerate and live in subspace of n basis vectors
    assertLoop(~res_ && size(subspace,2) == n,i);

    % same but using two inequalities instead of equality
    P_ = polytope(blkdiag(P.A,[1;-1]),[P.b;0;0]);
    [res_,subspace] = isFullDim(P_);
    assertLoop(~res_ && size(subspace,2) == n,i)

end

% ------------------------------ END OF CODE ------------------------------
