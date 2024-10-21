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

% empty zonotopes
Z = zonotope.empty(2);
assert(representsa(Z,'emptySet') && dim(Z) == 2);
Z = zonotope(zeros(3,0));
assert(representsa(Z,'emptySet'))
assert(size(Z.c,1) == 3);
assert(size(Z.G,1) == 3);
Z = zonotope(zeros(3,0),[]);
assert(representsa(Z,'emptySet'))
assert(size(Z.c,1) == 3);
assert(size(Z.G,1) == 3);
Z = zonotope(zeros(3,0),zeros(3,0));
assert(representsa(Z,'emptySet'))
assert(size(Z.c,1) == 3);
assert(size(Z.G,1) == 3);


% random center, random generator matrix
c = [3; 3; 2];
G = [2 -4 -6 3 5; 1 -7 3 -5 2; 0 4 -7 3 2];
Zmat = [c,G];

% admissible initializations
Z = zonotope(c,G);
assert(compareMatrices(Z.c,c) && compareMatrices(Z.G,G));

Z = zonotope(c);
assert(compareMatrices(Z.c,c) && isempty(Z.G) && size(G,1) == 3);

Z = zonotope(Zmat);
assert(compareMatrices(Z.c,Zmat(:,1)))
assert(compareMatrices(Z.G,Zmat(:,2:end)));


% wrong instantiations
c_plus1 = [4; 6; -2; 3];
G_plus1 = [2 -4 -6 3 5; 1 -7 3 -5 2; 0 4 -7 3 2; 2 0 5 -4 2];
randLogicals = randn(size(G)) > 0;
c_NaN = c; c_NaN(2) = NaN;
G_NaN = G; G_NaN(randLogicals) = NaN;

% center and generator matrix do not match
assertThrowsAs(@zonotope,'CORA:wrongInputInConstructor',c_plus1,G);
assertThrowsAs(@zonotope,'CORA:wrongInputInConstructor',c,G_plus1);

% center is empty
assertThrowsAs(@zonotope,'CORA:wrongInputInConstructor',[],G);

% center has NaN entry
assertThrowsAs(@zonotope,'CORA:wrongValue',c_NaN,G);

% generator matrix has NaN entries
assertThrowsAs(@zonotope,'CORA:wrongValue',c,G_NaN);

% too many input arguments
assertThrowsAs(@zonotope,'CORA:numInputArgsConstructor',c,G,G);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
