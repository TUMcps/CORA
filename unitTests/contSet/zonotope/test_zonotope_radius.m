function res = test_zonotope_radius
% test_zonotope_radius - unit test function of radius
%
% Syntax:  
%    res = test_zonotope_radius
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

% create a zonotope
Z_cent = zeros(2,1);
Z_gens = [1 3;
          2 1];
Z = zonotope([Z_cent, Z_gens]);
% plot(Z);

% compute radius
r = radius(Z);

% analytical solution
r_true = 5;

% check results
tol = 1e-9;
res = abs(r - r_true) < tol;

if res
    disp('test_radius successful');
else
    disp('test_radius failed');
end

%------------- END OF CODE --------------