function res = test_zonotope_norm
% test_zonotope_norm - unit test function of norm
%
% Syntax:  
%    res = test_zonotope_norm
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

% Author:       Mark Wetzlinger, Victor Gassmann
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

TOL = 1e-6;
res = true;

% instantiate zonotope
c = zeros(2,1);
G = [2 5 4 3; -4 -6 2 3];
Z = zonotope(c,G);


% 2-norm test
val2_exact = norm(Z,2,'exact');
val2_ub = norm(Z,2,'ub');
val2_ubc = norm(Z,2,'ub_convex');
V = vertices(Z);
if (val2_exact-val2_ub) > TOL || (val2_exact-val2_ubc) > TOL ...
        || abs(val2_exact-max(sqrt(sum(V.^2))))/val2_exact > TOL
    res = false;
end


if res
    disp('test_zonotope_norm successful');
else
    disp('test_zonotope_norm failed');
end

%------------- END OF CODE --------------
