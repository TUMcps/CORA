function res = test_interval_enclosePoints
% test_interval_enclosePoints - unit test function of enclosePoints
%
% Syntax:
%    res = test_interval_enclosePoints
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
% Written:       27-July-2021
% Last update:   03-December-2023 (MW, add unbounded case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% point cloud
pts = [1 4 2 5 3 2 4 3 2 5; ...
       6 8 9 7 6 9 8 6 8 7];
I = interval.enclosePoints(pts);
I_true = interval([1; 6],[5; 9]);
res(end+1,1) = isequal(I,I_true);

% unbounded
pts = [-Inf -2 1 5; ...
        2    3 4 -1; ...
        1    3 1 Inf];
I = interval.enclosePoints(pts);
I_true = interval([-Inf;-1;1],[5;4;Inf]);
res(end+1,1) = isequal(I,I_true);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
