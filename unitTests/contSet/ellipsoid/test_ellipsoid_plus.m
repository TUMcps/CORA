function res = test_ellipsoid_plus
% test_ellipsoid_plus - unit test function of plus
%
% Syntax:  
%    res = test_ellipsoid_plus
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
% Written:      14-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nPoints = 100;
nRuns = 2;
for i=10:5:15
    for j=1:nRuns
        E1 = ellipsoid.generateRandom(false,i);
        E2 = ellipsoid.generateRandom(false,i);
        E = E1+E2;
        Y1 = sample(E1,nPoints);
        Y2 = sample(E2,nPoints);
        %check if Y1+Y2 \in E
        for k=1:length(Y1)
            Y = repmat(Y1(:,k),1,nPoints)+Y2;
            if ~all(containsPoint(E,Y))
                res = false;
                break;
            end
        end
    end
end

if res
    disp('test_ellipsoid_plus successful');
else
    disp('test_ellipsoid_plus failed');
end
%------------- END OF CODE --------------
