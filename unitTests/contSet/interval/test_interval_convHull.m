function res = test_interval_convHull
% test_interval_convHull - unit test function of the convex hull
%
% Syntax:
%    res = test_interval_convHull
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
% Written:       12-March-2021
% Last update:   03-December-2023 (MW, add bounded and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case: convHull is empty set
I_empty = interval.empty(3);
I_3D = interval(-rand(3,1),rand(3,1));
res(end+1,1) = isequal(convHull(I_3D,I_empty),I_empty);

% bounded, 1D
I1 = interval(-2,4);
% ...I2 contained in I1
I2 = interval(-1,3);
I_convHull = convHull(I1,I2);
res(end+1,1) = isequal(I_convHull,I1);
% ...I2 intersects I1
I2 = interval(3,6);
I_convHull = convHull(I1,I2);
I_true = interval(-2,6);
res(end+1,1) = isequal(I_convHull,I_true);
% ...I2 is disjoint from I1
I2 = interval(5,6);
I_convHull = convHull(I1,I2);
I_true = interval(-2,6);
res(end+1,1) = isequal(I_convHull,I_true);

% unbounded, 2D
I1 = interval([-2;-Inf],[2;1]);
I2 = interval([-1;3],[Inf;8]);
I_convHull = convHull(I1,I2);
I_true = interval([-2;-Inf],[Inf;8]);
res(end+1,1) = isequal(I_convHull,I_true);


% dimension mismatch
I_dim3 = interval(-rand(3,1),rand(3,1));
I_dim7 = interval(-rand(7,1),rand(7,1));
try
    convHull(I_dim3,I_dim7);
    res = false;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
