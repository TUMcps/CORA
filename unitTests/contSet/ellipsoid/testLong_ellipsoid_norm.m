function res = testLong_ellipsoid_norm
% testLong_ellipsoid_norm - unit test function of norm
%
% Syntax:
%    res = testLong_ellipsoid_norm
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
nrOfTests = 1000;

for i=1:nrOfTests
    
    % random dimension
    n = randi(30);
    
    % ellipsoid is a ball: norm = 1
    Q = eye(n);
    q = zeros(n,1);
    E = ellipsoid(Q,q);
    
    % check result
    if ~withinTol(norm(E),1)
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
