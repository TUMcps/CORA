function res = test_conHyperplane_dim
% test_conHyperplane_dim - unit test function of dimension
%
% Syntax:
%    res = test_conHyperplane_dim
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

% 2D, empty
hyp = conHyperplane.empty(2);
res(end+1,1) = dim(hyp) == 2;

% 2D, normal vector and offset
a = [1 2]; b = 2;
hyp = conHyperplane(a,b);
res(end+1,1) = dim(hyp) == 2;

% 1D, normal vector, offset, and constraints
a = 1; b = 1; C = 1; d = 2;
hyp = conHyperplane(a,b,C,d);
res(end+1,1) = dim(hyp) == 1;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
