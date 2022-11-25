function res = testLongDuration_ellipsoid_enclose
% testLongDuration_ellipsoid_enclose - unit test function of enclose
%
% Syntax:  
%    res = testLongDuration_ellipsoid_enclose
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
% Written:      13-March-2019
% Last update:  22-March-2021 (extended to 1 degenerate ellipsoid)
% Last revision:---

%------------- BEGIN CODE --------------

%We need at least one non-degenerate ellipsoid to ensure
%the union of E1,E2 is not the empty set.
nRuns = 2;
res = true;
bools = [false,true];
for i=10:5:15
    for j=1:nRuns
        for k=1:2
            try
            E1 = ellipsoid.generateRandom(false,i);
            E2 = ellipsoid.generateRandom(bools(k),i);
            E = enclose(E1,E2);
            if ~in(E,E1) || ~in(E,E2)
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
if res
    disp('testLongDuration_ellipsoid_enclose successful');
else
    disp('testLongDuration_ellipsoid_enclose failed');
end


%------------- END OF CODE --------------
