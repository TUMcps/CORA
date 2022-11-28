function res = test_conZonotope_cubMap
% test_conZonotope_cubMap - unit test function for cubic multiplication of 
%                           constrained zonotopes
%
% Syntax:  
%    res = test_conZonotope_cubMap
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
% Written:      30-October-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

%% ANALYTICAL TESTS

% TEST 1: Mixed Multiplication

% define zonotope
Z = [0 1 -1; 1 2 0];
A = [1 1];
b = 0;
cZ = conZonotope(Z,A,b);

% define third-order tensor
temp = [1 -1; 0 2];
T{1,1} = temp;
T{1,2} = temp;
T{2,1} = temp;
T{2,2} = temp;

% compute cubic map
cZres = cubMap(cZ,cZ,cZ,T);

% define ground truth
temp = [2 3 1 4 7 6 1 -1 -2 1 9 3 -3 12 -1 -4 21 3 -3 -7 3 -1 1 -1];
Z_ = [temp;temp];
A_ = zeros(3,23);
A_(1,1) = 1;
A_(1,2) = 1;
A_(2,3) = 1;
A_(3,5) = 1;
A_(3,8) = 1;
b_ = [0;0;0];

% check for correctness
if ~compareMatrices(Z_,cZres.Z) || ...
        ~compareMatrices(A_,cZres.A) || ~all(withinTol(b_,cZres.b))
    res = false;
end



% TEST 2: Cubic Multiplication

% define zonotope
Z = [0 1 -1; 1 2 0];
A = [1 1];
b = 0;
cZ = conZonotope(Z,A,b);

% define third-order tensor
temp = [1 -1; 0 2];
T{1,1} = temp;
T{1,2} = temp;
T{2,1} = temp;
T{2,2} = temp;

% compute cubic map
cZres = cubMap(cZ,T);

% define ground truth
temp = [16 13 -1 14 -4 21 -7 3 -1];
Z_ = [temp;temp];
A_ = [1 1 0 0 0 0 0 0];
b_ = 0;

% check for correctness
if ~compareMatrices(Z_,cZres.Z) || ...
        ~compareMatrices(A_,cZres.A) || ~all(withinTol(b_,cZres.b))
    res = false;
end

%------------- END OF CODE --------------