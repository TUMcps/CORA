function res = test_polyZonotope_relaxExponents
% test_polyZonotope_relaxExponents - unit test function for exponent relaxation
%
% Syntax:
%    res = test_polyZonotope_relaxExponents
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

% Authors:       Tobias Ladner
% Written:       25-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test simple cases
pZ = polyZonotope([1;1]);
pZ_relax = pZ.relaxExponents(1);
assert(isequal(pZ,pZ_relax));

pZ = polyZonotope(zonotope([1;1],[2 4; 5 2]));
pZ_relax = pZ.relaxExponents(1);
assert(isequal(pZ,pZ_relax));

% test paper
pZ = polyZonotope([2;2],[1 -1; 2 3],[],[1 1; 4 2]);
pZ_relax = pZ.relaxExponents(1);
pZ_exp = polyZonotope([2;2],[0;5],[0.25;0.5],[1;2]);
assert(isequal(pZ_relax,pZ_exp));

% test graphs
pZ = polyZonotope([2;2],[1 -1; 2 3],[],[1 1; 4 2]);
[~,G] = pZ.relaxExponents(1,'greedy');
assert(isa(G,'digraph'));

pZ = polyZonotope([2;2],[1 -1; 2 3],[],[1 1; 4 2]);
[~,G] = pZ.relaxExponents(1,'all');
assert(isa(G,'graph'));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
