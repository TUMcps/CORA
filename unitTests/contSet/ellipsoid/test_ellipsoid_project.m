function res = test_ellipsoid_project
% test_ellipsoid_project - unit test function of project
%
% Syntax:
%    res = test_ellipsoid_project
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

% Authors:       Mark Wetzlinger
% Written:       27-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create a random psd Q matrix
n = 4;
O = orth(randn(n));
D = diag(abs(randn(n,1)) + 0.3);
Q = O*D*O';

% generate ellipsoid
E = ellipsoid(Q);

% project ellipsoid
projDim = [2 3];
E_proj1 = project(E,projDim);

% true solution
Q_proj = Q(projDim(1):projDim(2),projDim(1):projDim(2));
E_true = ellipsoid(Q_proj);

% logical indexing
projDim = [false true true];
E_proj2 = project(E,projDim);

% check properties
res(1) = all(all(withinTol(E_true.Q,E_proj1.Q)));
res(2) = E_true.dim == E_proj1.dim;
res(3) = all(all(withinTol(E_true.Q,E_proj2.Q)));
res(4) = E_true.dim == E_proj2.dim;
% summary of checks
res = all(res);

% ------------------------------ END OF CODE ------------------------------
