function res = testLongDuration_component_ellipsoid_distanceEllipsoid
% testLongDuration_component_ellipsoid_distanceEllipsoid - unit test 
% function of distanceEllipsoid
%
% Syntax:  
%    res = testLongDuration_component_ellipsoid_distanceEllipsoid
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

% Author:       Victor Gassmann
% Written:      18-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
bools = [false,true];
for i=5:5:10
    for j=1:nRuns
        for k=1:2 
            try
                E1 = ellipsoid.generateRandom(i,false);
                E2 = ellipsoid.generateRandom(i,bools(k));
                N = 2*i;
                m = ceil(N*rand);
                B = randPoint(E1,N,'extreme');
                b = B(:,m);
                s = E1.q + rand*(b-E1.q);
                % E1 and E2 are intersecting
                E2 = ellipsoid(E2.Q,s);
                if ~withinTol(distance(E1,E2),0,E1.TOL)
                    res = false;
                    break;
                end
                IntE1 = interval(E1);
                q3 = E1.q+2.1*rad(IntE1);
                [U,S,~] = svd(E1.Q);
                r = ceil(i*rand);
                S(r,r) = bools(k)*S(r,r);
                % guaranteed to not intersect with E1
                E3 = ellipsoid(U*S*U',q3);
                if withinTol(distance(E1,E3),0,E1.TOL)
                    res = false;
                    break;
                end
            catch ME
                if strcmp(ME.identifier,'CORA:solverIssue')
                    disp('Randomly generated ellipsoids caused solver issues! Ignoring...');
                    continue;
                end
                rethrow(ME);
            end
        end
        if ~res
            break;
        end
    end
    if ~res
        break;
    end
end
%------------- END OF CODE --------------