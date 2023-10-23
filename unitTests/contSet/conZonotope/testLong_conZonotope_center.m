function res = testLong_conZonotope_center
% testLong_conZonotope_center - unit test function for center computation
%
% Syntax:
%    res = testLong_conZonotope_center
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

% Authors:       Mark Wetzlinger
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nrTests = 100;

for i=1:nrTests
    % random dimension
    n = randi([2,8]);

    % random number of generators
    nrGens = randi(n*5);

    % random number of constraints
    nrCon = randi(n);

    % random center of zonotope
    c = randn(n,1);

    % instantiate random constrained zonotope
    cZ = conZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,...
        'NrConstraints',nrCon,'Center',c);

    % instantiate random constrained zonotope without constraints
    cZ_Z = conZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,...
        'NrConstraints',0,'Center',c);

    % compute centers
    c_cZ = center(cZ);
    c_Z = center(cZ_Z);

    % center of constrained zonotope without constraints has to be center
    % of zonotope
    if ~all(withinTol(c_Z,c,0))
        res = false; break
    end

    % check containment of center of actually constrained zonotope
    if ~contains(cZ,c_cZ)
        res = false; break
    end

end

% ------------------------------ END OF CODE ------------------------------
