function res = test_polyZonotope_reduce
% test_polyZonotope_reduce - unit test function of order reduction
%
% Syntax:
%    res = test_polyZonotope_reduce
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

% Authors:       Niklas Kochdumper
% Written:       29-March-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% TEST 1

% create polynomial zonotope
c = [0;0];
G = [0 2 1;3 0 1];
E = [1 0 1;0 1 1];
GI = [0;0];
pZ = polyZonotope(c,G,GI,E);

% reduce the polynomial zonotope
pZred = reduce(pZ,'girard',1);

% define ground truth
GI =  [3 0; 0 4];
c = [0;0];

% check for correctness
assert(all(withinTol(c,pZred.c)))
assert(compareMatrices(GI,pZred.GI))
assert(isempty(pZred.G))
assert(isempty(pZred.E))
assert(isempty(pZred.id))


% TEST 2

% create polynomial zonotope
c = [0;0];
G = [0 3 1 1; 2 0 1 -3];
E = [1 0 1 2;0 1 1 0];
GI = [-1 -4; -2 -1];
pZ = polyZonotope(c,G,GI,E);

% reduce the polynomial zonotope
pZred = reduce(pZ,'girard',2);

% define ground truth
c = [0.5;-1.5];
G = [3;0];
GI = [-4 2.5 0; -1 0 6.5];
E = 1;

% check for correctness
assert(all(withinTol(c,pZred.c)))
assert(compareMatrices(G,pZred.G))
assert(compareMatrices(GI,pZred.GI))

end

% ------------------------------ END OF CODE ------------------------------
