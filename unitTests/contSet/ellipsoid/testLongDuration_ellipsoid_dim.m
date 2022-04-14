function res = testLongDuration_ellipsoid_dim
% testLongDuration_ellipsoid_dim - unit test function of dim
%
% Syntax:  
%    res = testLongDuration_ellipsoid_dim
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

% empty case: dim = 0
res_empty = true;
E = ellipsoid();
if dim(E) ~= 0
    res_empty = false;
end


% random tests
res_rand = true;
nrOfTests = 100;

for i=1:nrOfTests
    
    % random dimension
    n = randi(15);
    
    % non-degenerate case
    E = ellipsoid.generateRandom(n,false);
    % check result
    if dim(E) ~= n
        res_rand = false; break;
    end
    
    % degenerate case
    E = ellipsoid.generateRandom(n,true);
    % check result
    if dim(E) ~= n
        res_rand = false; break;
    end

end

% combine results
res = res_empty && res_rand;

if res
    disp('testLongDuration_ellipsoid_dim successful');
else
    disp('testLongDuration_ellipsoid_dim failed');
end

%------------- END OF CODE --------------