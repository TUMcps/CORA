function res = test_ellipsoid_interval
% test_ellipsoid_interval - unit test function of interval
%
% Syntax:  
%    res = test_ellipsoid_interval
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

% Author:       Victor Gaﬂmann
% Written:      16-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
for i=10:5:15
    for j=1:nRuns
        E = ellipsoid.generateRandom(false,i);
        Y = sample(E,1000);
        %compute interval
        I = interval(E);
        %check if all points are in the interval
        if ~all(containsPoint(I,Y))
            res = false;
            break;
        end
    end
    if ~res
        break;
    end
end

if res
    disp('test_ellipsoid_interval successful');
else
    disp('test_ellipsoid_interval failed');
end
%------------- END OF CODE --------------
