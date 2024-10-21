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

% empty case
I = interval.empty(2);
assert(contains(I,I));

% zonotope-in-interval
I = interval([-3;-2],[5;4]);
Z_in = zonotope([0.5;0], [2, 1; 1,-0.7]);
assert(contains(I, Z_in));
Z_out = zonotope([6.5;-3],[2, 1; 1,-0.7]);
assert(~contains(I, Z_out));

% point-containment case
I = interval([-3; -9; -4; -7; -1], [4;   2;  6;  3;  8]);
p_inside = [0, 0, 0, 0, 0;
            1,-4, 3,-6, 5;
           -2,-6,-2, 2, 7;
           -3, 2, 6,-7, 8]';
res = contains(I,p_inside);
assert(numel(res) == 4)
assert(all(res));
p_outside = [5, 3, 7, 4, 9;
             1,-4,-5,-6, 5;
            -2, 3,-2, 6, 7;
            -3, 2, 6,-7,10]';
assert(all(~contains(I,p_outside)));

% unbounded intervals
I1 = interval(-Inf,Inf);
I2 = interval(-Inf,0);
assert(contains(I1,I2));
I2 = interval(0,Inf);
assert(contains(I1,I2));
p = [-Inf,0,Inf];
assert(all(contains(I1,p)));

% n-d arrays
lb = [];
lb(:,:,1,1) = [1 2; 3 5];
lb(:,:,1,2) = [0 -1; -2 3];
lb(:,:,1,3) = [1 1; -1 0];
lb(:,:,2,1) = [-3 2; 0 1];
ub = [];
ub(:,:,1,1) = [1.5 4; 4 10];
ub(:,:,1,2) = [1 2; 0 4];
ub(:,:,1,3) = [2 3; -0.5 2];
ub(:,:,2,1) = [-1 3; 0 2];
I = interval(lb,ub);
c = center(I);
res = contains(I,cat(5,c,c));
assert(isequal(size(res),[1,2]));
assert(all(res))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
