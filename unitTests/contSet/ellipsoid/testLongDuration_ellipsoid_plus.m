function res = testLongDuration_ellipsoid_plus
% testLongDuration_ellipsoid_plus - unit test function of plus
%
% Syntax:  
%    res = testLongDuration_ellipsoid_plus
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
% Written:      14-October-2019
% Last update:  16-March-2021
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
rng(123456);
for i=10:5:15
    for j=1:nRuns
        try
            E1 = ellipsoid.generateRandom(i,false);
            E2 = ellipsoid.generateRandom(i,true);
            E3 = ellipsoid.generateRandom(i); 
            % compute outer approx (check overloaded syntax)
            Eo1 = E1+E2;
            % check regular syntax
            Eo2 = plus(E1,{E2,E3},'o');
            % compute inner approx
            Ei1 = plus(E1,E2,'i');
            Ei2 = plus(E1,{E2,E3},'i');
            % check if inner contained in outer
            if ~in(Eo1,Ei1) || ~in(Eo2,Ei2)
                res = false;
                break;
            end
            Y1 = randPoint(E1,i,'extreme');
            Y2 = randPoint(E2,i,'extreme');
            Y3 = randPoint(E3,i,'extreme');
            %check if Yi+Yj \in E
            Yres1 = sumPoints(Y1,Y2);
            Yres2 = sumPoints(Y1,Y2,Y3);
            if ~in(Eo1,Yres1) || ~in(Eo2,Yres2)
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

if res
    disp('testLongDuration_ellipsoid_plus successful');
else
    disp('testLongDuration_ellipsoid_plus failed');
end
%------------- END OF CODE --------------