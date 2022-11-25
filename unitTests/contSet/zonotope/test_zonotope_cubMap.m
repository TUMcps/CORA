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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Niklas Kochdumper
% Written:      16-August-2018
% Last update:  01-May-2020 (MW, cubicMultiplication -> cubMap)
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% ANALYTICAL TESTS

% TEST 1: Mixed Multiplication

% define zonotope
Z = [0 1 -1; 1 2 0];
zono = zonotope(Z);

% define third-order tensor
temp = [1 -1; 0 2];
T{1,1} = temp;
T{1,2} = temp;
T{2,1} = temp;
T{2,2} = temp;

% compute cubic map
Zres = cubMap(zono,zono,zono,T);

% define ground truth
temp = [2 3 1 4 7 1 0 -1 1 6 9 3 12 21 3 0 -3 3 -2 -3 -1 -4 -7 -1 0 1 -1];
Z_ = [temp;temp];

% check for correctness
if any(any(Z_-Zres.Z))
    error('zonotope/cubMap: analytical test (mixed mul.) failed!');
end



% TEST 2: Cubic Multiplication

% define zonotope
Z = [0 1 -1; 1 2 0];
zono = zonotope(Z);

% define third-order tensor
temp = [1 -1; 0 2];
T{1,1} = temp;
T{1,2} = temp;
T{2,1} = temp;
T{2,2} = temp;

% compute cubic map
Zres = cubMap(zono,T);

% define ground truth
temp = [16 13 -1 14 -4 0 21 -7 3 -1];
Z_ = [temp;temp];

% check for correctness
if any(any(Z_-Zres.Z))
    error('zonotope/cubMap: analytical test failed!');
end




res = true;

end

%------------- END OF CODE --------------