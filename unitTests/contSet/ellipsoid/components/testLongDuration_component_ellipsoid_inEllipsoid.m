function res = testLongDuration_component_ellipsoid_inEllipsoid
% testLongDuration_component_ellipsoid_inEllipsoid - unit test function of inEllipsoid
%
% Syntax:  
%    res = testLongDuration_component_ellipsoid_inEllipsoid
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
% Written:      16-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
for i=10:5:15
    for j=1:nRuns
        try
            E1 = ellipsoid.generateRandom(i,true);
            E2 = ellipsoid.generateRandom(i,false);
            % check whether E1 dg, E2 non-d results in E2\subseteq E1 = false
            if in(E1,E2)
                res = false;
                break;
            end
            % E3 contains E2
            E3 = ellipsoid(0.9*E2.Q,E2.q);
            % check if E3 contains E2
            if ~in(E2,E3)
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
%------------- END OF CODE --------------