function res = test_interval_contains
% test_interval_contains - unit test function of contains
%
% Syntax:
%    res = test_interval_contains
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

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       27-September-2019
% Last update:   12-March-2021 (MW, add empty case)
%                01-July-2021 (AK, integrated test_interval_containsPoint)
%                03-December-2023 (MW, add unbounded case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
I = interval.empty(2);
res(end+1,1) = contains(I,I);

% zonotope-in-interval
I = interval([-3;-2],[5;4]);
Z_in = zonotope([0.5;0], [2, 1; 1,-0.7]);
res(end+1,1) = contains(I, Z_in);
Z_out = zonotope([6.5;-3],[2, 1; 1,-0.7]);
res(end+1,1) = ~contains(I, Z_out);

% point-containment case
I = interval([-3; -9; -4; -7; -1], [4;   2;  6;  3;  8]);
p_inside = [0, 0, 0, 0, 0;
            1,-4, 3,-6, 5;
           -2,-6,-2, 2, 7;
           -3, 2, 6,-7, 8]';
res(end+1,1) = all(contains(I,p_inside));
p_outside = [5, 3, 7, 4, 9;
             1,-4,-5,-6, 5;
            -2, 3,-2, 6, 7;
            -3, 2, 6,-7,10]';
res(end+1,1) = all(~contains(I,p_outside));

% unbounded intervals
I1 = interval(-Inf,Inf);
I2 = interval(-Inf,0);
res(end+1,1) = contains(I1,I2);
I2 = interval(0,Inf);
res(end+1,1) = contains(I1,I2);
p = [-Inf,0,Inf];
res(end+1,1) = all(contains(I1,p));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
