function res = test_polyZonotope_zonotope
% test_polyZonotope_zonotope - unit test function for zonotope
%    over-approximation of a polynomial zonotope
%
% Syntax:  
%    res = test_polyZonotope_zonotope
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
% Written:      26-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% ANALYTICAL TESTS

% TEST 1

% create polynomial zonotope
c = [1;2];
G = [1 2 1 -3; 1 -1 2 -1];
expMat = [1 0 0 2; 0 1 2 1];
Grest = [];
pZ = polyZonotope(c,G,Grest,expMat);

% reduce the polynomial zonotope
zono = zonotope(pZ);
c_ = zono.Z(:,1);
G_ = zono.Z(:,2:end);

% define ground truth
c = [1.5; 3];
G =  [1 2 0.5 -3; 1 -1 1 -1];

% check for correctness
if any(c-c_)
    error('test_polyZonotope_zonotope: analytical test 1 failed!');
end
for i = 1:size(G,2)
   if ~ismember(G(:,i)',G_)
       error('test_polyZonotope_zonotope: analytical test 1 failed!');
   end
end





res = true;

%------------- END OF CODE --------------