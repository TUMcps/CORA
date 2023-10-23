function res = test_nn_nnHelper_calcSquared()
% test_nn_nnHelper_calcSquared - tests the nnHelper.calcSquared function
%
% Syntax:
%    res = test_nn_nnHelper_calcSquared()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnHelper/calcSquared, nnHelper/calcSquaredE

% Authors:       Tobias Ladner
% Written:       17-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
c = [1];
G = [2 1.5 1];
GI = []; % only using GI=[] to avoid over-approximation
E = [1 0 3;0 1 1];

% compute pZ^2
[c2, G2, GI2] = nnHelper.calcSquared( ...
    c, G, GI, E, ...
    c, G, GI, E, true);
E2 = nnHelper.calcSquaredE(E, E, true);
E2 = [E, E, E2];
pZ2 = polyZonotope(c2, G2, GI2, E2);

% compute pZ^3
[c3, G3, GI3] = nnHelper.calcSquared( ...
    c2, G2, GI2, E2, ...
    c, G, GI, E, false);
E3 = nnHelper.calcSquaredE(E2, E, false);
E3 = [E2, E, E3];
pZ3 = polyZonotope(c3, G3, GI3, E3);

% compute reference
pZ_ref = polyZonotope(c,G,GI,E);
pZ2_ref = quadMap(pZ_ref, {1});
pZ3_ref = cubMap(pZ_ref, {1});

% check equality
res = isequal(pZ2, pZ2_ref);
res = res && isequal(pZ3, pZ3_ref);

end

% ------------------------------ END OF CODE ------------------------------
