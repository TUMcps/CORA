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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      27-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create a random psd Q matrix
dim = 4;
O = orth(randn(dim));
D = diag(abs(randn(dim,1)) + 0.3);
Q = O*D*O';

% generate ellipsoid
E = ellipsoid(Q);

% project ellipsoid
projDim = [2 3];
E_proj = project(E,projDim);

% true solution
Q_proj = Q(projDim(1):projDim(2),projDim(1):projDim(2));
E_true = ellipsoid(Q_proj);

% check properties
tol = 1e-9;
res(1) = all(all(abs(E_true.Q - E_proj.Q) < tol));
res(2) = E_true.dim == E_proj.dim;
res(3) = E_true.isdegenerate == E_proj.isdegenerate;
% summary of checks
res = all(res);


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end

%------------- END OF CODE --------------