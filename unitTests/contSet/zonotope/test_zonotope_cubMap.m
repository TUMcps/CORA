function res = test_zonotope_cubMap
% test_zonotope_cubMap - unit test function for cubic multiplication of 
%                        zonotopes
%
% Syntax:
%    res = test_zonotope_cubMap
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

% Authors:       Niklas Kochdumper
% Written:       16-August-2018
% Last update:   01-May-2020 (MW, cubicMultiplication -> cubMap)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% TEST 1: Mixed Multiplication

% define zonotope
Z = [0 1 -1; 1 2 0];
zono = zonotope(Z);

% define third-order tensor
tensorSlice = [1 -1; 0 2];
T{1,1} = tensorSlice;
T{1,2} = tensorSlice;
T{2,1} = tensorSlice;
T{2,2} = tensorSlice;

% compute cubic map
Zres = cubMap(zono,zono,zono,T);

% define ground truth
truthRow = [2 3 1 4 7 1 0 -1 1 6 9 3 12 21 3 0 -3 3 -2 -3 -1 -4 -7 -1 0 1 -1];
Z_ = [truthRow;truthRow];

% check for correctness
assert(compareMatrices(Z_,[Zres.c, Zres.G]))


% TEST 2: Cubic Multiplication

% define zonotope
Z = [0 1 -1; 1 2 0];
zono = zonotope(Z);

% define third-order tensor
tensorSlice = [1 -1; 0 2];
T{1,1} = tensorSlice;
T{1,2} = tensorSlice;
T{2,1} = tensorSlice;
T{2,2} = tensorSlice;

% compute cubic map
Zres = cubMap(zono,T);

% define ground truth
truthRow = [16 13 -1 14 -4 0 21 -7 3 -1];
Z_ = [truthRow;truthRow];

% check for correctness
assert(compareMatrices(Z_,[Zres.c,Zres.G]))

end

% ------------------------------ END OF CODE ------------------------------
