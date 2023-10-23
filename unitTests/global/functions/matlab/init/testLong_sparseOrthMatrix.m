function res = testLong_sparseOrthMatrix
% testLong_sparseOrthMatrix - unit test function for instantiation of
%    sparse orthogonal matrices
%
% Syntax:
%    res = testLong_sparseOrthMatrix()
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true, wait for failure
res = true;

nrTests = 1000;
for i=1:nrTests
    % random dimension
    n = randi(20);

    % init sparse orthogonal matrix
    Q = sparseOrthMatrix(n);

    % check orthogonality
    if ~all(withinTol(vecnorm(Q),1))
        res = false; break
    end

end

% ------------------------------ END OF CODE ------------------------------
