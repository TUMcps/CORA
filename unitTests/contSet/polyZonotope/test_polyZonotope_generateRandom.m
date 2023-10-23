function res = test_polyZonotope_generateRandom
% test_polyZonotope_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_polyZonotope_generateRandom
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
% Written:       19-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty call
pZ = polyZonotope.generateRandom();

% values for tests
n = 3;
nrGens = 10;
nrFac = 3;
nrIndGens = 5;

% only dimension
pZ = polyZonotope.generateRandom('Dimension',n);
res = dim(pZ) == n;

% only number of generators (no check)
pZ = polyZonotope.generateRandom('NrGenerators',nrGens);

% only number of factors
pZ = polyZonotope.generateRandom('NrFactors',nrFac);
res(end+1,1) = size(pZ.E,1) == nrFac;

% only number of independent generators
pZ = polyZonotope.generateRandom('NrIndGenerators',nrIndGens);
res(end+1,1) = size(pZ.GI,2) == nrIndGens;

% dimension and number of generators
pZ = polyZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens);
res(end+1,1) = dim(pZ) == n;

% dimension and number of factors
pZ = polyZonotope.generateRandom('Dimension',n,'NrFactors',nrFac);
res(end+1,1) = dim(pZ) == n && size(pZ.E,1) == nrFac;

% dimension, number of generators, and number of factors
pZ = polyZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,'NrFactors',nrFac);
res(end+1,1) = dim(pZ) == n && size(pZ.E,1) == nrFac;

% dimension, number of factors, and number of independent generators
pZ = polyZonotope.generateRandom('Dimension',n,'NrIndGenerators',nrIndGens,'NrFactors',nrFac);
res(end+1,1) = dim(pZ) == n && size(pZ.E,1) == nrFac;

% dimension, number of generators, and number of independent generators
pZ = polyZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,'NrIndGenerators',nrIndGens);
res(end+1,1) = dim(pZ) == n && size(pZ.GI,2) == nrIndGens;

% dimension, number of generators, number of factors, and number of independent generators
pZ = polyZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,...
    'NrFactors',nrFac,'NrIndGenerators',nrIndGens);
res(end+1,1) = dim(pZ) == n && size(pZ.E,1) == nrFac && size(pZ.GI,2) == nrIndGens;

% unify results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
