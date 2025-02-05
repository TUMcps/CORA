function res = test_polyZonotope_polyMap
% test_polyZonotope_polyMap - unit test function of polyMap
%
% Syntax:
%    res = test_polyZonotope_polyMap
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
% Written:       23-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% test 1

% polynomial zonotope
pZ = polyZonotope(0,1,[],[2;3]);

% compute polynomial map
pZ_res = polyMap(pZ,1,2);

% check correctness
assert(pZ_res.c == 0)
assert(pZ_res.G == 1)
assert(all(pZ_res.E == [4;6],"all"))

%% test 2

% polynomial zonotope
pZ = polyZonotope([1;2],[1 -2 1; 2 3 1],[0;0],[1 0 2;0 1 1]);

% polynomial map [x(2)^3 + 2*x(1)*x(2)^2; -x(2)^3]
coeff = [1 2;-1 0];
E = [0 1;3 2];

% compute polynomial map
pZ_res = polyMap(pZ,coeff,E);

c = [ 16 ; -8 ];
G = [ 48 44 48 88 24 16 72 24 -9 56 56 28 56 21 16 16 17 3 ; -24 -36 -24 -72 -54 -8 -48 -54 -27 -24 -36 -12 -36 -27 -6 -6 -9 -1 ];
GI = zeros(2,0);
E = [ 1 0 2 1 0 3 2 1 0 3 2 4 3 2 4 5 4 6 ; 0 1 0 1 2 0 1 2 3 1 2 1 2 3 2 2 3 3 ];
id = [ 1 ; 2 ];
pZ_true = polyZonotope(c,G,GI,E,id);

assert(isequal(pZ_res,pZ_true))

%% test 3 - all-zero exponents

coeff = [1 4 2 2;-1 -3 0 1];
E = [0 0 1 0;3 0 2 0];

% compute polynomial map
pZ_res = polyMap(pZ,coeff,E);

c = [ 22 ; -10 ];
G = [ 48 44 48 88 24 16 72 24 -9 56 56 28 56 21 16 16 17 3 ; -24 -36 -24 -72 -54 -8 -48 -54 -27 -24 -36 -12 -36 -27 -6 -6 -9 -1 ];
GI = zeros(2,0);
E = [ 1 0 2 1 0 3 2 1 0 3 2 4 3 2 4 5 4 6 ; 0 1 0 1 2 0 1 2 3 1 2 1 2 3 2 2 3 3 ];
id = [ 1 ; 2 ];
pZ_true = polyZonotope(c,G,GI,E,id);

assert(isequal(pZ_res,pZ_true))


%% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
