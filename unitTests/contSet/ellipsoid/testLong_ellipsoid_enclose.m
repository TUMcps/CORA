function res = testLong_ellipsoid_enclose
% testLong_ellipsoid_enclose - unit test function of enclose
%
% Syntax:
%    res = testLong_ellipsoid_enclose
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
% Written:       13-March-2019
% Last update:   22-March-2021 (extended to 1 degenerate ellipsoid)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%We need at least one non-degenerate ellipsoid to ensure
%the union of E1,E2 is not the empty set.
nRuns = 20;
res = true;
bools = [false,true];
for i=10:5:15
    for j=1:nRuns
        for k=1:2
            try
            %%% generate all variables necessary to replicate results
            E1 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
            E2 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            %%%
            E = enclose(E1,E2);
            if ~contains(E,E1) || ~contains(E,E2)
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
end

% ------------------------------ END OF CODE ------------------------------
