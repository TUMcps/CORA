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

res = true;

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
if ~all(withinTol(c,c_)) || ~compareMatrices(G,G_)
    res = false;
end

%------------- END OF CODE --------------