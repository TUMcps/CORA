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
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      30-July-2018
% Last update:  06-December-2022 (TL: more tests)
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% TEST 1

% create polynomial zonotope
c = [3;4];
G = [2 0 1;0 2 1];
expMat = [1 0 1;0 1 1];
Grest = [0;0];
pZ = polyZonotope(c,G,Grest,expMat);

% calculate enclosing interval
I = interval(pZ,'bnb');

% define ground truth
I_ = interval([0;1],[6;7]);

% check for correctness
if ~isequal(I,I_)
    res = false; return
end

% test empty set
pZ_e = polyZonotope();
if supportFunc(pZ_e,[1;1],'upper') ~= -Inf ...
        || supportFunc(pZ_e,[1;1],'lower') ~= +Inf
    res = false; return
end

% test [lower, upper] = range
dir = [1 1];
range = supportFunc(pZ, dir, 'range');
lower = supportFunc(pZ, dir, 'lower');
upper = supportFunc(pZ, dir, 'upper');

if ~isequal(range, interval(lower, upper))
    res = false; return
end


% TEST 2 (check point containment)

c = [0;-1];
G = [-0.5 0.4 0.04 -0.04; 0.6 1.9 1.7 -0.2];
Grest = [];
E = [0 4 2 1; 1 0 3 2];
pZ = polyZonotope(c, G, Grest, E);

% slightly under-approximative
ground_truth = interval([-0.5059;-2.1168], [0.9; 3.34]);

methods = ["interval", 'split', 'bnb', 'bnbAdv'];
for method = methods
    I = cartProd( ...
      supportFunc(pZ, [1 0], 'range', method), ...
      supportFunc(pZ, [0 1], 'range', method) ...
    );

    if ~I.contains(ground_truth, 'exact', 1e-15)
        res = false; return
    end
end

%------------- END OF CODE --------------