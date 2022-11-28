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
for i=[2,10:5:15]
    for j=1:nRuns
        try
            %%% generate all variables necessary to replicate results
            E1 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
            E2 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',true);
            E3 = ellipsoid.generateRandom('Dimension',i); 
            Y1 = randPoint(E1,i,'extreme');
            Y2 = randPoint(E2,i,'extreme');
            Y3 = randPoint(E3,i,'extreme');
            %%%
            % compute outer approx (check overloaded syntax)
            Eo1 = E1+E2;
            % check regular syntax
            Eo2 = plus(E1,[E2,E3],'outer');
            % 'inner' only exists for 'L' syntax
%             % compute inner approx
%             Ei1 = plus(E1,E2,'inner');
%             Ei2 = plus(E1,[E2,E3],'inner');
%             % check if inner contained in outer
%             if ~contains(Eo1,Ei1) || ~contains(Eo2,Ei2)
%                 res = false;
%                 break;
%             end
            %check if Yi+Yj \in E
            Yres1 = sumPoints(Y1,Y2);
            Yres2 = sumPoints(Y1,Y2,Y3);
            if ~contains(Eo1,Yres1) || ~contains(Eo2,Yres2)
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
    path = pathFailedTests(mfilename());
    save(path,'E1','E2','E3','Y1','Y2','Y3');
end

%------------- END OF CODE --------------