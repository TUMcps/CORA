function res = testLong_ellipsoid_generateRandom
% testLong_ellipsoid_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = testLong_ellipsoid_generateRandom
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
% Written:       13-March-2021
% Last update:   19-March-2021 (VG, removed false check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% random tests
res = true;
nrOfTests = 100;

for i=1:nrOfTests
    
    % random dimension
    n = randi(15);
    E = ellipsoid.generateRandom('IsDegenerate',false);
    assertLoop(isFullDim(E),i)

    % ellipsoid with random dimension, Q, and center, degenerate
    E = ellipsoid.generateRandom('IsDegenerate',true);
    assertLoop(~isFullDim(E),i)
    
    % fixed dimension, non-degenerate
    E = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',false);
    assertLoop(dim(E) == n,i)
    assertLoop(isFullDim(E),i)

    % fixed dimension, degenerate
    E = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',true);
    assertLoop(dim(E) == n,i)
    assertLoop(~isFullDim(E),i)

end

% ------------------------------ END OF CODE ------------------------------
