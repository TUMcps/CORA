function res = testLong_ellipsoid_contains
% testLong_ellipsoid_contains - unit test function of contains
%
% Syntax:
%    res = testLong_ellipsoid_contains
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
% See also: none

% Authors:       Victor Gassmann
% Written:       15-October-2019
% Last update:   07-August-2020
%                19-March-2021
% Last revision: 22-September-2024 (MW, integrate components)

% ------------------------------ BEGIN CODE -------------------------------

% containment of ellipsoids
assert(aux_containsEllipsoid());

% containment of zonotopes
assert(aux_containsZonotope());

% point containment
assert(aux_containsPoint());


% combine results
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_containsEllipsoid()

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
            assertLoop(~contains(E1,E2),i,j)

            % E3 contains E2
            E3 = ellipsoid(0.5*E2.Q,(1+1e-4)*E2.q);
            % check if E3 contains E2
            assertLoop(contains(E2,E3),i,j)

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

function res = aux_containsZonotope()

res = true;
nRuns = 2;
for i=2:3
    for j=1:nRuns
        % do not test vertices in high dims as comp. complexity prevents it
        % for high dimensions
        Z = zonotope.generateRandom('Dimension',2);
        E_Zo = ellipsoid(Z,'outer:exact');
        assertLoop(contains(E_Zo,Z),i,j)

        E_Zi = ellipsoid(Z,'inner:exact');
        assertLoop(contains(E_Zo,E_Zi),i,j)
    end
end

end

function res = aux_containsPoint()

res = true;
runs = 10;
bools = [false,true];
for i=1:5:30
    for j=1:runs
        for k=1:2
            %%% generate all variables necessary to replicate results
            N = 10*i;
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            samples = randPoint(E,N);
            %%%
            
            assertLoop(all(contains(E,samples)),i,j)
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
