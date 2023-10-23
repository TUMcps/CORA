function res = testLong_ellipsoid_volume
% testLong_ellipsoid_volume - unit test function of volume
%
% Syntax:
%    res = testLong_ellipsoid_volume
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
    %%% generate all variables necessary to replicate results
    % random dimension
    n = randi(15);
    % degenerate case: ball of dim n in space (n+1)
    E_d = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',true);
    %%%

    % ellipsoid is a ball
    Q = eye(n);
    q = zeros(n,1);
    E = ellipsoid(Q,q);
    vol_ball = (pi^(n/2) / gamma(1+n/2)) * 1^n;
    % check result
    if ~withinTol(vol_ball,volume(E))
        res = false; break;
    end
    
    % check result
    try
        volume(E_d);
    catch
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
