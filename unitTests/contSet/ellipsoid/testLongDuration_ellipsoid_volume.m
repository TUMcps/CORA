function res = testLongDuration_ellipsoid_volume
% testLongDuration_ellipsoid_volume - unit test function of volume
%
% Syntax:  
%    res = testLongDuration_ellipsoid_volume
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

% Author:       Mark Wetzlinger
% Written:      13-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;

% empty case: volume?
res_empty = true;
E = ellipsoid();
if abs(volume(E)) > tol
    res_empty = false;
end


% random tests
res_rand = true;
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
        res_rand = false; break;
    end
    
    % check result
    try
        volume(E_d);
    catch
        res_rand = false; break;
    end

end

% combine results
res = res_empty && res_rand;

if ~res
    path = pathFailedTests(mfilename());
    save(path,'n','E_d');
end

%------------- END OF CODE --------------