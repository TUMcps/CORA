function res = testLong_spectraShadow_generateRandom
% testLong_spectraShadow_generateRandom - unit test function of 
%    spectraShadow.generateRandom
%
% Syntax:
%    res = testLong_spectraShadow_generateRandom
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

% Authors:       Adrian Kulmburg
% Written:       05-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty call
SpS = spectraShadow.generateRandom();

% values for tests
n = 2;
m = 5;
c = [2;1];

% only dimension
SpS = spectraShadow.generateRandom('Dimension',n,'NrGenerators',m);
assert((dim(SpS) == n) & (size(SpS.G,2) == m));

% only center
SpS = spectraShadow.generateRandom('Center',c,'NrGenerators',m);
assert(all(abs(SpS.c - c) < eps));

% only degeneracy = true
SpS = spectraShadow.generateRandom('IsDegenerate',true,'NrGenerators',m);
assert(~isFullDim(SpS));

% only degeneracy = false
SpS = spectraShadow.generateRandom('IsDegenerate',false,'NrGenerators',m);
assert(isFullDim(SpS));

% dimension and center
SpS = spectraShadow.generateRandom('Dimension',n,'Center',c,'NrGenerators',m);
assert(dim(SpS) == n && all(abs(SpS.c - c) < eps));

% dimension and degeneracy = true
SpS = spectraShadow.generateRandom('Dimension',n,'IsDegenerate',true,'NrGenerators',m);
assert(dim(SpS) == n && ~isFullDim(SpS));

% dimension and degeneracy = false
SpS = spectraShadow.generateRandom('Dimension',n,'IsDegenerate',false,'NrGenerators',m);
assert(dim(SpS) == n && isFullDim(SpS));

% center and degeneracy = true
SpS = spectraShadow.generateRandom('Center',c,'IsDegenerate',true,'NrGenerators',m);
assert(all(abs(SpS.c - c) < eps) && ~isFullDim(SpS));

% center and degeneracy = false
SpS = spectraShadow.generateRandom('Center',c,'IsDegenerate',false,'NrGenerators',m);
assert(all(abs(SpS.c - c) < eps) && isFullDim(SpS));

% dimension, center, and degeneracy = true
SpS = spectraShadow.generateRandom('Dimension',n,'Center',c,'IsDegenerate',true,'NrGenerators',m);
assert(dim(SpS) == n && all(abs(SpS.c - c) < eps) && ~isFullDim(SpS));

% dimension, center, and degeneracy = false
SpS = spectraShadow.generateRandom('Dimension',n,'Center',c,'IsDegenerate',false,'NrGenerators',m);
assert(dim(SpS) == n && all(abs(SpS.c - c) < eps) && isFullDim(SpS));


% unify results
res = true;

% ------------------------------ END OF CODE ------------------------------
