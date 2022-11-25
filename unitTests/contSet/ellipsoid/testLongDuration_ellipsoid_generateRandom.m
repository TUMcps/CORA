function res = testLongDuration_ellipsoid_generateRandom
% testLongDuration_ellipsoid_generateRandom - unit test function of generateRandom
%
% Syntax:  
%    res = testLongDuration_ellipsoid_generateRandom
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
% Last update:  19-March-2021 (VG: removed false check)
% Last revision:---

%------------- BEGIN CODE --------------

% random tests
res = true;
nrOfTests = 100;

for i=1:nrOfTests
    
    E = ellipsoid.generateRandom(false);
    if ~isFullDim(E)
        res = false; break;
    end
    % ellipsoid with random dimension, Q, and center, degenerate
    E = ellipsoid.generateRandom(true);
    if isFullDim(E)
        res = false; break;
    end
    
    % random dimension
    n = randi(15);
    
    % fixed dimension, non-degenerate
    E = ellipsoid.generateRandom(n,false);
    if dim(E) ~= n || ~isFullDim(E)
        res = false; break;
    end
    % fixed dimension, degenerate
    E = ellipsoid.generateRandom(n,true);
    if dim(E) ~= n || isFullDim(E)
        res = false; break;
    end

end


if res
    disp('testLongDuration_ellipsoid_generateRandom successful');
else
    disp('testLongDuration_ellipsoid_generateRandom failed');
end

%------------- END OF CODE --------------