function res = testLong_ellipsoid_plus
% testLong_ellipsoid_plus - unit test function of plus
%
% Syntax:
%    res = testLong_ellipsoid_plus
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       14-October-2019
% Last update:   16-March-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nRuns = 2;
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
            Eo = E1+E2;
            % 'inner' only exists for 'L' syntax
%             % compute inner approx
%             Ei = plus(E1,E2,'inner');
%             % check if inner contained in outer
%             if ~contains(Eo,Ei)
%                 res = false;
%                 break;
%             end
            %check if Yi+Yj \in E
            Yres = sumPoints(Y1,Y2);
            if ~contains(Eo,Yres)
                res = false; return
            end
        catch ME
            if strcmp(ME.identifier,'CORA:solverIssue')
                disp('Randomly generated ellipsoids caused solver issues! Ignoring...');
                continue;
            end
            rethrow(ME);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
