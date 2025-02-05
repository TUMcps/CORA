function res = test_zonotope_randPoint
% test_zonotope_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_zonotope_randPoint
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

% Authors:       Mark Wetzlinger
% Written:       05-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-7;

% nD, empty zonotope
n = 4;
Z = zonotope.empty(n);
p = randPoint(Z);
assert(all(size(p) == [n,0]));


% 2D, zonotope that's a single point
Z = zonotope([2;3]);
p = randPoint(zonotope([2;3]),1);
p_extr = randPoint(zonotope([2;3]),1,'extreme');
assert(size(p,2) == 1 && compareMatrices(p,Z.c,tol));
assert(size(p_extr,2) == 1 && compareMatrices(p,Z.c,tol));


% 3D, degenerate zonotope
Z = zonotope([1;2;-1],[1 3 -2; 1 0 1; 2 3 -1]);
numPoints = 10;
p_random = randPoint(Z,numPoints,'standard');
assert(all(contains(Z,p_random,'exact',tol)));
p_random = randPoint(Z,numPoints,'uniform');
assert(all(contains(Z,p_random,'exact',tol)));
p_random = randPoint(Z,numPoints,'uniform:hitAndRun');
assert(all(contains(Z,p_random,'exact',tol)));
p_random = randPoint(Z,numPoints,'uniform:billiardWalk');
assert(all(contains(Z,p_random,'exact',tol)));
% less points than extreme points
nrPtsExtreme = ceil(2^n * 0.5);
p_random = randPoint(Z,nrPtsExtreme,'extreme');
assert(all(contains(Z,p_random,'exact',tol)));
% as many points as extreme points
nrPtsExtreme = ceil(2^n);
p_random = randPoint(Z,nrPtsExtreme,'extreme');
assert(all(contains(Z,p_random,'exact',tol)));
% more points than extreme points
nrPtsExtreme = ceil(2^n * 5);
p_random = randPoint(Z,nrPtsExtreme,'extreme');
assert(all(contains(Z,p_random,'exact',tol)));
% Gaussian sampling (not supported for degenerate...)
% p_random = randPoint(Z,10,'gaussian',0.8);


% 3D, parallelotope
Z = zonotope([1;2;-1],[1 3 -2; 1 0 1; 4 1 -1]);
numPoints = 10;
p_random = randPoint(Z,numPoints,'standard');
assert(all(contains(Z,p_random,'exact',tol)));
p_random = randPoint(Z,numPoints,'uniform');
assert(all(contains(Z,p_random,'exact',tol)));
p_random = randPoint(Z,numPoints,'uniform:hitAndRun');
assert(all(contains(Z,p_random,'exact',tol)));
p_random = randPoint(Z,numPoints,'uniform:billiardWalk');
assert(all(contains(Z,p_random,'exact',tol)));
% less points than extreme points
nrPtsExtreme = ceil(2^n * 0.5);
p_random = randPoint(Z,nrPtsExtreme,'extreme');
assert(all(contains(Z,p_random,'exact',tol)));
% as many points as extreme points
nrPtsExtreme = ceil(2^n);
p_random = randPoint(Z,nrPtsExtreme,'extreme');
assert(all(contains(Z,p_random,'exact',tol)));
% more points than extreme points
nrPtsExtreme = ceil(2^n * 5);
p_random = randPoint(Z,nrPtsExtreme,'extreme');
assert(all(contains(Z,p_random,'exact',tol)));
% Gaussian sampling (containment not guaranteed)
p_random = randPoint(Z,10,'gaussian',0.8);


% 3D, general zonotope
Z = zonotope([0;1;1],[0 -1 1 4 2 3 -2; 0 1 4 2 1 3 -8; 0 1 -3 2 1 2 6]);
numPoints = 10;
p_random = randPoint(Z,numPoints,'standard');
assert(all(contains(Z,p_random,'exact',tol)));
p_random = randPoint(Z,numPoints,'uniform');
assert(all(contains(Z,p_random,'exact',tol)));
p_random = randPoint(Z,numPoints,'uniform:hitAndRun');
assert(all(contains(Z,p_random,'exact',tol)));
p_random = randPoint(Z,numPoints,'uniform:billiardWalk');
assert(all(contains(Z,p_random,'exact',tol)));
% less points than extreme points
nrPtsExtreme = ceil(2^n * 0.5);
p_random = randPoint(Z,nrPtsExtreme,'extreme');
assert(all(contains(Z,p_random,'exact',tol)));
% as many points as extreme points
nrPtsExtreme = ceil(2^n);
p_random = randPoint(Z,nrPtsExtreme,'extreme');
assert(all(contains(Z,p_random,'exact',tol)));
% more points than extreme points
nrPtsExtreme = ceil(2^n * 5);
p_random = randPoint(Z,nrPtsExtreme,'extreme');
assert(all(contains(Z,p_random,'exact',tol)));
% Gaussian sampling (containment not guaranteed)
p_random = randPoint(Z,10,'gaussian',0.8);


% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
