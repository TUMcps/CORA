function res = testLongDuration_ellipsoid_norm
% testLongDuration_ellipsoid_norm - unit test function of norm
%
% Syntax:  
%    res = testLongDuration_ellipsoid_norm
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

% empty case: norm = -Inf
res_empty = true;
E = ellipsoid();
if norm(E) ~= -Inf
    res_empty = false;
end


% random tests
res_rand = true;
nrOfTests = 1000;

for i=1:nrOfTests
    
    % random dimension
    n = randi(30);
    
    % ellipsoid is a ball: norm = 1
    Q = eye(n);
    q = zeros(n,1);
    E = ellipsoid(Q,q);
    
    % check result
    if abs(1 - norm(E)) > tol
        res_rand = false; break;
    end

end

% combine results
res = res_empty && res_rand;

if ~res
    path = pathFailedTests(mfilename());
    save(path,'n','E');
end

%------------- END OF CODE --------------