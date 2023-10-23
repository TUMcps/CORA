function res = test_ellipsoid_generateRandom
% test_ellipsoid_generateRandom - unit test function of 
%    ellipsoid.generateRandom
%
% Syntax:
%    res = test_ellipsoid_generateRandom
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

% Authors:       Victor Gassmann, Mark Wetzlinger
% Written:       26-July-2021
% Last update:   19-May-2022 (name-value pair syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty call
E = ellipsoid.generateRandom();

% values for tests
n = 5;
q = [2;1;-1;3;4];

% only dimension
E = ellipsoid.generateRandom('Dimension',n);
res = dim(E) == n;

% only center
E = ellipsoid.generateRandom('Center',q);
res(end+1,1) = all(abs(E.q - q) < eps);

% only degeneracy = true
E = ellipsoid.generateRandom('IsDegenerate',true);
res(end+1,1) = ~isFullDim(E);

% only degeneracy = false
E = ellipsoid.generateRandom('IsDegenerate',false);
res(end+1,1) = isFullDim(E);

% dimension and center
E = ellipsoid.generateRandom('Dimension',n,'Center',q);
res(end+1,1) = dim(E) == n && all(abs(E.q - q) < eps);

% dimension and degeneracy = true
E = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',true);
res(end+1,1) = dim(E) == n && ~isFullDim(E);

% dimension and degeneracy = false
E = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',false);
res(end+1,1) = dim(E) == n && isFullDim(E);

% center and degeneracy = true
E = ellipsoid.generateRandom('Center',q,'IsDegenerate',true);
res(end+1,1) = all(abs(E.q - q) < eps) && ~isFullDim(E);

% center and degeneracy = false
E = ellipsoid.generateRandom('Center',q,'IsDegenerate',false);
res(end+1,1) = all(abs(E.q - q) < eps) && isFullDim(E);

% dimension, center, and degeneracy = true
E = ellipsoid.generateRandom('Dimension',n,'Center',q,'IsDegenerate',true);
res(end+1,1) = dim(E) == n && all(abs(E.q - q) < eps) && ~isFullDim(E);

% dimension, center, and degeneracy = false
E = ellipsoid.generateRandom('Dimension',n,'Center',q,'IsDegenerate',false);
res(end+1,1) = dim(E) == n && all(abs(E.q - q) < eps) && isFullDim(E);


% unify results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
