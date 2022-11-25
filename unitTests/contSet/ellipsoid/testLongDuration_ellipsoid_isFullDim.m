function res = testLongDuration_ellipsoid_isFullDim
% testLongDuration_ellipsoid_isFullDim - unit test function of isFullDim
%
% Syntax:  
%    res = testLongDuration_ellipsoid_isFullDim
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

% empty case: not full-dimensional
res_empty = true;
E = ellipsoid();
if isFullDim(E)
    res_empty = false;
end


% random tests
res_rand = true;
nrOfTests = 100;

for i=1:nrOfTests
    
    % random dimension
    n = randi(15);
    
    % non-degenerate case: full-dimensional
    E = ellipsoid.generateRandom(n,false);
    % check result
    if ~isFullDim(E)
        res_rand = false; break;
    end
    
    % degenerate case: not full-dimensional
    E = ellipsoid.generateRandom(n,true);
    % check result
    if isFullDim(E)
        res_rand = false; break;
    end

end

% combine results
res = res_empty && res_rand;

if res
    disp('testLongDuration_ellipsoid_isFullDim successful');
else
    disp('testLongDuration_ellipsoid_isFullDim failed');
end

%------------- END OF CODE --------------