function res = testLong_ellipsoid_dim
% testLong_ellipsoid_dim - unit test function of dim
%
% Syntax:
%    res = testLong_ellipsoid_dim
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% random tests
res = true;
nrOfTests = 100;

for i=1:nrOfTests
    
    % random dimension
    n = randi(15);
    
    % non-degenerate case
    E = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',false);
    % check result
    if dim(E) ~= n
        res = false; break;
    end
    
    % degenerate case
    E = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',true);
    % check result
    if dim(E) ~= n
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
