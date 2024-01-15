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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty zonotopes
Z = zonotope.empty(2);
res(end+1,1) = representsa(Z,'emptySet') && dim(Z) == 2;
Z = zonotope(zeros(3,0));
res(end+1,1) = representsa(Z,'emptySet') ...
    && size(Z.c,1) == 3 && size(Z.G,1) == 3;
Z = zonotope(zeros(3,0),[]);
res(end+1,1) = representsa(Z,'emptySet') ...
    && size(Z.c,1) == 3 && size(Z.G,1) == 3;
Z = zonotope(zeros(3,0),zeros(3,0));
res(end+1,1) = representsa(Z,'emptySet') ...
    && size(Z.c,1) == 3 && size(Z.G,1) == 3;


% random center, random generator matrix
c = [3; 3; 2];
G = [2 -4 -6 3 5; 1 -7 3 -5 2; 0 4 -7 3 2];
Zmat = [c,G];

% admissible initializations
Z = zonotope(c,G);
res(end+1,1) = compareMatrices(Z.c,c) && compareMatrices(Z.G,G);

Z = zonotope(c);
res(end+1,1) = compareMatrices(Z.c,c) && isempty(Z.G) && size(G,1) == 3;

Z = zonotope(Zmat);
res(end+1,1) = compareMatrices(Z.c,Zmat(:,1)) ...
    && compareMatrices(Z.G,Zmat(:,2:end));


% combine results
res = all(res);

% wrong initializations
if CHECKS_ENABLED

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

end

% ------------------------------ END OF CODE ------------------------------
