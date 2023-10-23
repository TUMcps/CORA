function res = test_nn_nnHelper_computeBoundsPolyZono()
% test_nn_nnHelper_computeBoundsPolyZono - tests the 
%     nnHelper.compBoundsPolyZono function
%
% Syntax:
%    res = test_nn_nnHelper_computeBoundsPolyZono()
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
% See also: nnHelper/compBoundsPolyZono

% Authors:       Tobias Ladner
% Written:       17-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
c = [0];
G = [2 0 1];
GI = [0.5 0.1];
E = [1 0 3;0 1 1];
pZ = polyZonotope(c, G, GI, E);

ind = find(prod(ones(size(E))-mod(E, 2), 1) == 1);
ind_ = setdiff(1:size(E, 2), ind);

% approximate bounds
[l,u] = nnHelper.compBoundsPolyZono(c, G, GI, E, ind, ind_, true);
bounds = interval(l, u);
bounds_ref = interval(zonotope(pZ));

res = isequal(bounds, bounds_ref);

% tighter bounds using splitting
[l,u] = nnHelper.compBoundsPolyZono(c, G, GI, E, ind, ind_, false);
bounds = interval(l, u);
bounds_ref = interval(pZ, 'split');

res = res && isequal(bounds, bounds_ref);

end

% ------------------------------ END OF CODE ------------------------------
