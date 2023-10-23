function res = testLong_component_ellipsoid_inEllipsoid
% testLong_component_ellipsoid_inEllipsoid - unit test function of
%    inEllipsoid
%
% Syntax:
%    res = testLong_component_ellipsoid_inEllipsoid
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
% Written:       16-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nRuns = 2;
for i=10:5:15
    for j=1:nRuns
        try
            %%% generate all variables necessary to replicate results
            E1 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',true);
            E2 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
            %%%
            
            % check whether E1 dg, E2 non-d results in E2\subseteq E1 = false
            if contains(E1,E2)
                res = false;
                return;
            end
            % E3 contains E2
            E3 = ellipsoid(0.5*E2.Q,(1+1e-4)*E2.q);
            % check if E3 contains E2
            if ~contains(E2,E3)
                res = false;
                return;
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
