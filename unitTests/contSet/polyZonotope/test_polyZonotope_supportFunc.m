function res = test_polyZonotope_supportFunc
% test_polyZonotope_supportFunc - unit test function for calculating bounds
%    of the polynomial zonotope along a specific direction 
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

res = true;

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

% empty set
pZ_e = polyZonotope();
if supportFunc(pZ_e,[1;1],'upper') ~= -Inf || supportFunc(pZ_e,[1;1],'lower') ~= +Inf
    res = false;
end

% check for correctness
if ~isequal(inter,inter_)
    res = false;
end


%------------- END OF CODE --------------