function res = test_interval_radius
% test_interval_radius - unit test function of radius
%
% Syntax:
%    res = test_interval_radius
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
% Written:       27-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create interval (all positive)
Int = interval([-2;-3], [4;1]);

% compute radius
Int_rad = radius(Int);

% true solution
true_rad = sqrt(13);

res = Int_rad == true_rad;

% ------------------------------ END OF CODE ------------------------------
