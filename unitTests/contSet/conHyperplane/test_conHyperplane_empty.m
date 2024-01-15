function res = test_conHyperplane_empty
% test_conHyperplane_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_conHyperplane_empty
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
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D
hyp = conHyperplane.empty(2);
res(end+1,1) = dim(hyp) == 2;
res(end+1,1) = representsa(hyp,'emptySet');

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
