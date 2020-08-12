function res = test_ellipsoid_and
% test_ellipsoid_and - unit test function of plus
%
% Syntax:  
%    res = test_ellipsoid_and
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
for i=10:5:20
    for j=1:nRuns
        E1 = ellipsoid.generateRandom(false,i);
        E2 = ellipsoid.generateRandom(false,i);
        E = E1&E2;
        Y1 = sample(E1,nPoints);
        %see which parts of Y1 are both E1 and E2
        res1 = containsPoint(E1,Y1);
        res2 = containsPoint(E2,Y1);
        %find points which are in the intersection
        R = res1 & res2;
        if isempty(E)
            if any(R)
                res = false;
                break;
            end
        else
            Y = Y1(:,R);
            %check if Y points are in E
            if ~all(containsPoint(E,Y))
                res = false;
                break;
            end
        end
    end
end

if res
    disp('test_ellipsoid_and successful');
else
    disp('test_ellipsoid_and failed');
end
%------------- END OF CODE --------------
