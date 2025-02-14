function res = test_zonotope_generateRandom
% test_zonotope_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_zonotope_generateRandom
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
% Written:       27-September-2019
% Last update:   19-May-2022 (name-value pair syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty call
Z = zonotope.generateRandom();

% values for tests
n = 3;
c = [2;1;-1];
nrGens = 10;
type = 'exp';

% only dimension
Z = zonotope.generateRandom('Dimension',n);
assert(dim(Z) == n);

% only center
Z = zonotope.generateRandom('Center',c);
assert(compareMatrices(c,Z.c));

% only number of generators
Z = zonotope.generateRandom('NrGenerators',nrGens);
assert(size(Z.G,2) == nrGens);

% only type (no check)
Z = zonotope.generateRandom('Distribution',type);

% dimension and number of generators
Z = zonotope.generateRandom('Dimension',n,'NrGenerators',nrGens);
assert(dim(Z) == n && size(Z.G,2) == nrGens);

% center and number of generators
Z = zonotope.generateRandom('Center',c,'NrGenerators',nrGens);
assert(compareMatrices(c,Z.c) && size(Z.G,2) == nrGens);

% center and type
Z = zonotope.generateRandom('Center',c,'Distribution',type);
assert(compareMatrices(c,Z.c));

% number of generators and type
Z = zonotope.generateRandom('NrGenerators',nrGens,'Distribution',type);
assert(size(Z.G,2) == nrGens);

% center, number of generators, and type
Z = zonotope.generateRandom('Center',c,'NrGenerators',nrGens,'Distribution',type);
assert(compareMatrices(c,Z.c) && size(Z.G,2) == nrGens);


% unify results
res = true;

% ------------------------------ END OF CODE ------------------------------
