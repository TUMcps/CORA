function res = testLong_ellipsoid_isFullDim
% testLong_ellipsoid_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = testLong_ellipsoid_isFullDim
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

% empty case: not full-dimensional
res = true;
% E = ellipsoid.empty(2);
% if isFullDim(E)
%     res = false;
% end


% random tests
nrOfTests = 100;

for i=1:nrOfTests
    
    % random dimension
    n = randi(15);
    
    % non-degenerate case: full-dimensional
    E = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',false);
    % check result
    if ~isFullDim(E)
        res = false; break;
    end
    
    % degenerate case: not full-dimensional
    E = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',true);
    % check result
    if isFullDim(E)
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
