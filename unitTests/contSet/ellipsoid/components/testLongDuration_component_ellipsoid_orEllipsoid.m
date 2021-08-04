function res = testLongDuration_component_ellipsoid_orEllipsoid
% testLongDuration_component_ellipsoid_orEllipsoid - unit test function of 
% orEllipsoidIA
%
% Syntax:  
%    res = testLongDuration_component_ellipsoid_orEllipsoid
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
% Written:      17-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
bools = [false,true];
for i=10:5:15
    for j=1:nRuns
        for k=1:2
            try
                E1 = ellipsoid.generateRandom(i,false);
                E2 = ellipsoid.generateRandom(i,false);
                E3 = ellipsoid.generateRandom(i,bools(k)); 
                % compute outer approx (check overloaded syntax)
                Eo1 = E1|E2;
                % check regular syntax
                Eo2 = or(E1,{E2,E3},'o');
                % compute inner approx
                Ei1 = or(E1,E2,'i');
                Ei2 = or(E1,{E2,E3},'i');
                % check if inner contained in outer
                if ~in(Eo1,Ei1) || ~in(Eo2,Ei2)
                    res = false;
                    break;
                end
                Y1 = randPoint(E1,i,'extreme');
                Y2 = randPoint(E2,i,'extreme');
                Y3 = randPoint(E3,i,'extreme');
                if ~in(Eo1,Y1) || ~in(Eo1,Y2) || ~in(Eo2,Y1) ||...
                   ~in(Eo2,Y2) || ~in(Eo2,Y3)
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