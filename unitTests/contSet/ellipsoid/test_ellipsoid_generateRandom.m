function res = test_ellipsoid_generateRandom
% test_ellipsoid_generateRandom - unit test function of 
%    ellipsoid.generateRandom
%
% Syntax:  
%    res = test_ellipsoid_generateRandom
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

% Author:       Victor Gassmann, Mark Wetzlinger
% Written:      26-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
    
E = ellipsoid.generateRandom(false);
if ~isFullDim(E)
    res = false; 
end
% ellipsoid with random dimension, Q, and center, degenerate
E = ellipsoid.generateRandom(true);
if isFullDim(E)
    res = false; 
end

% random dimension
n = randi(15);

% fixed dimension, non-degenerate
E = ellipsoid.generateRandom(n,false);
if dim(E) ~= n || ~isFullDim(E)
    res = false;
end
% fixed dimension, degenerate
E = ellipsoid.generateRandom(n,true);
if dim(E) ~= n || isFullDim(E)
    res = false;
end


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end
%------------- END OF CODE --------------
