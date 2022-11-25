function res = test_polyZonotope_supportFunc
% test_polyZonotope_supportFunc - unit test function for calculating bounds
%    of the polynomial zontope along a specific direction 
%
% Syntax:  
%    res = test_polyZonotope_supportFunc
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
% Written:      30-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% ANALYTICAL TESTS

% TEST 1

% create polynomial zonotope
c = [3;4];
G = [2 0 1;0 2 1];
expMat = [1 0 1;0 1 1];
Grest = [0;0];
pZ = polyZonotope(c,G,Grest,expMat);

% calculate enclosing interval
inter = interval(pZ,'bnb');

% define ground truth
inter_ = interval([0;1],[6;7]);

% check for correctness
if any(supremum(inter) ~= supremum(inter_)) || any(infimum(inter) ~= infimum(inter_))
    error('test_polyZonotope_supportFunc: analytical test 1 failed!');
end



res = true;

%------------- END OF CODE --------------