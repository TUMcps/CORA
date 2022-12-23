function res = test_zonotope_zonotope
% test_zonotope_zonotope - unit test function of zonotope (constructor)
%
% Syntax:  
%    res = test_zonotope_zonotope
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

% Author:       Mark Wetzlinger
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-12;

% empty zonotope
Z = zonotope();
res = isempty(Z);


% random center, random generator matrix
c = [3; 3; 2];
G = [2 -4 -6 3 5; 1 -7 3 -5 2; 0 4 -7 3 2];
Zmat = [c,G];

% admissible initializations
Z = zonotope(c,G);
if ~compareMatrices(center(Z),c) || ~compareMatrices(generators(Z),G)
    res = false;
end

Z = zonotope(c);
if ~compareMatrices(center(Z),c) || ~isempty(generators(Z))
    res = false;
end

Z = zonotope(Zmat);
if ~compareMatrices(center(Z),Zmat(:,1)) ...
        || ~compareMatrices(generators(Z),Zmat(:,2:end))
    res = false;
end

% wrong initializations
c_plus1 = [4; 6; -2; 3];
G_plus1 = [2 -4 -6 3 5; 1 -7 3 -5 2; 0 4 -7 3 2; 2 0 5 -4 2];
randLogicals = randn(size(G)) > 0;
c_NaN = c; c_NaN(2) = NaN;
G_NaN = G; G_NaN(randLogicals) = NaN;

% center and generator matrix do not match
try
    Z = zonotope(c_plus1,G); % <- should throw error here
    res = false;
end
try
    Z = zonotope(c,G_plus1); % <- should throw error here
    res = false;
end

% center is empty
try
    Z = zonotope([],G); % <- should throw error here
    res = false;
end


% center has NaN entry
try
    Z = zonotope(c_NaN,G); % <- should throw error here
    res = false;
end

% generator matrix has NaN entries
try
    Z = zonotope(c,G_NaN); % <- should throw error here
    res = false;
end

% too many input arguments
try
    Z = zonotope(c,G,G); % <- should throw error here
    res = false;
end 

%------------- END OF CODE --------------