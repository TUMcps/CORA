function res = test_interval_zonotope
% test_interval_zonotope - unit test function of zonotope
%
% Syntax:  
%    res = test_interval_zonotope
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
% Written:      28-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% TEST 1: Analytical ------------------------------------------------------
% create interval
lowerLimits = [-2; -5];
upperLimits = [3; 1];
I = interval(lowerLimits, upperLimits);

% convert to zonotope
Z = zonotope(I);

% true zonotope
cent = [0.5; -2];
gens = [2.5 0.0;
        0.0 3.0];
Z_true = zonotope([cent, gens]);

% compare results
tol = 1e-9;
res_analytical = all(abs(center(Z) - center(Z_true)) < tol) && ...
    all(all(abs(generators(Z) - generators(Z_true)) < tol));

% -------------------------------------------------------------------------

% TEST 2: random ----------------------------------------------------------
% create random interval
dim = floor(2 + 8*rand(1));
lowerLimits = -3+3*rand(dim,1);
upperLimits = 3*rand(dim,1);
Irand = interval(lowerLimits, upperLimits);

% convert to zonotope
Zrand = zonotope(Irand);

% true zonotope
centrand_true = center(Irand);
gensrand_true = diag(rad(Irand));
Zrand_true = zonotope([centrand_true, gensrand_true]);

% compare results
res_rand = all(abs(center(Zrand) - center(Zrand_true)) < tol) && ...
    all(all(abs(generators(Zrand) - generators(Zrand_true)) < tol));
% -------------------------------------------------------------------------


% add results
res = res_analytical && res_rand;

if res
    disp('test_zonotope successful');
else
    disp('test_zonotope failed');
end

%------------- END OF CODE --------------