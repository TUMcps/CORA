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
%    res - boolean 
%
% Example: 
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
    
    % random dimension
    n = randi(15);
    
    % ellipsoid is a ball
    Q = eye(n);
    q = zeros(n,1);
    E = ellipsoid(Q,q);
    vol_ball = (pi^(n/2) / gamma(1+n/2)) * 1^n;
    % check result
    if abs(vol_ball - volume(E)) > tol
        res_rand = false; break;
    end
    
    % degenerate case: ball of dim n in space (n+1)
    E = ellipsoid.generateRandom(n,true);
    % check result
    try
        volume(E);
    catch
        res_rand = false; break;
    end

end

% combine results
res = res_empty && res_rand;

if res
    disp('testLongDuration_ellipsoid_volume successful');
else
    disp('testLongDuration_ellipsoid_volume failed');
end

%------------- END OF CODE --------------