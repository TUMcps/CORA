function res = test_interval_volume
% test_interval_volume - unit test function of volume
%
% Syntax:  
%    res = test_interval_volume
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
lowerLimits = [-2; -4; -3];
upperLimits = [3; 1; 2];
int = interval(lowerLimits, upperLimits);

% compute volume
vol = volume(int);

% true volume
vol_true = 125;

% compare results
tol = 1e-9;
res_analytical = abs(vol - vol_true) < tol;
% -------------------------------------------------------------------------

% TEST 2: random ----------------------------------------------------------
% create random interval
dim = floor(2 + 8*rand(1));
lowerLimits = -3+3*rand(dim,1);
upperLimits = 3*rand(dim,1);
intRand = interval(lowerLimits, upperLimits);

% compute volume
vol = volume(intRand);

% true volume
vol_true = prod(upperLimits - lowerLimits);

% compare results
res_rand = abs(vol - vol_true) < tol;
% -------------------------------------------------------------------------


% add results
res = res_analytical && res_rand;

if res
    disp('test_volume successful');
else
    disp('test_volume failed');
end

%------------- END OF CODE --------------