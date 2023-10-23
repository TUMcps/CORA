function res = testLong_polytope_box
% testLong_polytope_box - unit test function of box
%
% Syntax:
%    res = testLong_polytope_box
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
% Written:       03-December-2022
% Last update:   27-July-2023 (MW, more tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 25;

for i=1:nrTests

    % random dimension
    n = randi(10);

    % instantiate polytope
    P = polytope.generateRandom('Dimension',n,'nrConstraints',2*n);
    
    % compute box
    B = box(P);

    % box has to contain the polytope
    if ~contains(B,P,'exact',1e-14)
        res = false; break
    end

    % instantiate from vertices (one vertex per quadrant to ensure that
    % none are redundant)


end

% ------------------------------ END OF CODE ------------------------------
