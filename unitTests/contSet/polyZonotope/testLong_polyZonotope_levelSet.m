function res = testLong_polyZonotope_levelSet
% testLong_polyZonotope_levelSet - unit test function for level-set
%       enclosure of polynomial zonotopes
%
% Syntax:
%    res = testLong_polyZonotope_levelSet
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
% See also: polyZonotope/levelSet

% Authors:       Niklas Kochdumper
% Written:       30-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    tol = 0.01;

    % generate random polynomial zonotope
    n = randi([2,3]); nrGens = randi([1,5]); nrFac = randi([1,5]);
    nrIndGens = randi([0,3]);

    pZ = polyZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,...
            'NrFactors',nrFac,'nrIndGenerators',nrIndGens);

    % enclose by level set
    ls = levelSet(pZ);

    % check for correctness
    p = randPoint(pZ,1000);

    assert(all(contains(ls,p,'exact',tol)));

    % gather results
    res = true;

% ------------------------------ END OF CODE ------------------------------
