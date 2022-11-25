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
% Written:      27-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create interval (all positive)
Int = interval([-2;-3], [4;1]);

% compute radius
Int_rad = radius(Int);

% true solution
true_rad = sqrt(13);

res = Int_rad == true_rad;

if res
    disp('test_radius successful');
else
    disp('test_radius failed');
end

%------------- END OF CODE --------------