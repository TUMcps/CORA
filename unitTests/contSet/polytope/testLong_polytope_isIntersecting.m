function res = testLong_polytope_isIntersecting
% testLong_polytope_isIntersecting - unit test function of intersection
%    check
%
% Syntax:
%    res = testLong_polytope_isIntersecting
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
% Written:       27-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% number of tests
nrTests = 20;

for i=1:nrTests

    % random dimension
    n = randi(10);

    % init polytopes
    P1 = polytope.generateRandom('Dimension',n,'NrConstraints',n+2);
    P2 = polytope.generateRandom('Dimension',n,'NrConstraints',n+2);

    % shift both sets by center
    P1 = P1 - center(P1);
    P2 = P2 - center(P2);

    % both sets contain the origin
    if ~isIntersecting(P1,P2)
        throw(CORAerror('CORA:testFailed'));
    end

    % shift P1 by random point inside P1
    P1 = P1 - randPoint(P1);
    % rotate by random matrix
    [M,~,~] = svd(randn(n));
    P1_ = M * P1;

    % both sets contain the origin
    if ~isIntersecting(P1,P1_)
        throw(CORAerror('CORA:testFailed'));
    end

end

% ------------------------------ END OF CODE ------------------------------
