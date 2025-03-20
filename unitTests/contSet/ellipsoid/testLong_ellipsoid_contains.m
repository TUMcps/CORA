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
% See also: -
%

% Authors:       Victor Gassmann, Adrian Kulmburg
% Written:       15-October-2019
% Last update:   07-August-2020
%                19-March-2021
%                10-June-2023 (AK, added zono-in-ellipsoid tests)
%                20-January-2025 (AK, added full containment checks)
% Last revision: 22-September-2024 (MW, integrate components)

% ------------------------------ BEGIN CODE -------------------------------

% containment of ellipsoids
assert(aux_containsEllipsoid());

% containment of zonotopes
assert(aux_containsZonotope());

% point containment
assert(aux_containsPoint());

% Set up general set containment check
S = 2*ellipsoid(eye(2));
Sdeg = [1 0; 0 0] * S;
Sempty = ellipsoid.empty(2);

implementedSets = {'capsule','conPolyZono','conZonotope','interval','polytope',...
                    'zonoBundle','zonotope','ellipsoid','taylm',...
                    'polyZonotope','conPolyZono','spectraShadow'};
setsNonExact = {'taylm','conPolyZono','polyZonotope',...
                'spectraShadow'};

additionalAlgorithms = {};
additionalAlgorithms_specificSets = {};

% check containments
checkAllContainments(S, Sdeg, Sempty, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);


% combine results
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_containsEllipsoid()

tol = 1e-5;

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
            assertLoop(~contains(E1,E2,'exact',tol),i,j)

            % E3 contains E2
            E3 = ellipsoid(0.5*E2.Q,(1+1e-4)*E2.q);
            % check if E3 contains E2
            assertLoop(contains(E2,E3,'exact',tol),i,j)

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
if false
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

% Check zono-in-ellipsoid methods
runs = 1:3;
dims = 2:5;
gens = 5:10;
bools = [false,true];
for r=runs
    for d = dims
        for g = gens
            for k=1:2
                E = ellipsoid.generateRandom('Dimension',d,'IsDegenerate',bools(k));
                Z = zonotope.generateRandom('Dimension',d,'NrGenerators',g);
                
                
                % Computing the exact containment problem; different
                % methods of calling the function
                res_exact_solo = E.contains(Z); % Using the exact method
                [res_exact,cert_exact,scaling_exact] = E.contains(Z);
                if res_exact_solo ~= res_exact
                    res = false; return
                end
                
                % Computing the approx containment problem
                res_approx_solo = E.contains(Z,'approx');
                [res_approx,cert_approx,scaling_approx] = E.contains(Z,'approx');
                if res_approx_solo ~= res_approx
                    res = false; return
                end
                
                if scaling_exact > scaling_approx+1e-5 || scaling_approx > scaling_exact * sqrt(pi/2)+1e-5
                    res = false; return
                end
            end
        end
    end
end


end

% ------------------------------ END OF CODE ------------------------------
