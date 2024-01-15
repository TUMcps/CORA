function res = test_interval_gridPoints
% test_interval_gridPoints - unit test function of gridPoints
%
% Syntax:
%    res = test_interval_gridPoints
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
% Written:       08-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-12;

% to plot, use:
% figure; hold on; box on;
% plot(I,[1,2],'LineWidth',2);
% scatter(vals(1,:),vals(2,:),16,'r','filled');

% 1. empty case
n = 2;
I = interval.empty(n);
vals = gridPoints(I,5);
res_empty = isempty(vals) && isnumeric(vals) && size(vals,1) == n;


% 2. interval
% full-dimensional non-empty case: segments scalar, greater than 1
I_fullDim = interval([-1;-2],[1;4]);
segments = 3;
vals = gridPoints(I_fullDim,segments);
vals_true = [-1 -1 -1  0  0  0  1  1  1
             -2  1  4 -2  1  4 -2  1  4];
res_int = all(all(abs(vals - vals_true) < tol));

% full-dimensional non-empty case: segments vector, some are 1
segments = [3;1];
vals = gridPoints(I_fullDim,segments);
vals_true = [-1  0  1
              1  1  1];
res_int(end+1) = all(all(abs(vals - vals_true) < tol));
         
% full-dimensional non-empty case: segments vector, greater than 1
segments = [3;2];
vals = gridPoints(I_fullDim,segments);
vals_true = [-1 -1  0  0  1  1
             -2  4 -2  4 -2  4];
res_int(end+1) = all(all(abs(vals - vals_true) < tol));

% full-dimensional non-empty case: transposed, segments vector, greater than 1
segments = [3 2];
vals = gridPoints(I_fullDim',segments);
vals_true = [-1 -1  0  0  1  1
             -2  4 -2  4 -2  4]';
res_int(end+1) = all(all(abs(vals - vals_true) < tol));

% non-empty case, partially radius 0: segments scalar, greater than 1
I_rpart0 = interval([-1;-2],[1;-2]);
segments = 3;
vals = gridPoints(I_rpart0,segments);
vals_true = [-1  0  1
             -2 -2 -2];
res_int(end+1) = all(all(abs(vals - vals_true) < tol));

% non-empty case, partially radius 0: segments vector, greater than 1
segments = [3;2];
vals = gridPoints(I_rpart0,segments);
vals_true = [-1  0  1
             -2 -2 -2];
res_int(end+1) = all(all(abs(vals - vals_true) < tol));

% non-empty case, partially radius 0: transposed, segments vector, greater than 1
segments = [3 2];
vals = gridPoints(I_rpart0',segments);
vals_true = [-1 -2; 0 -2; 1 -2];
res_int(end+1) = all(all(abs(vals - vals_true) < tol));

% non-empty case, partially radius 0: segments vector, some are 1
I_rpart0_3D = interval([-1;-2;-4],[1;-2;1]);
segments = [3;2;1];
vals = gridPoints(I_rpart0_3D,segments);
vals_true = [-1.0000         0    1.0000
             -2.0000   -2.0000   -2.0000
             -1.5000   -1.5000   -1.5000];
res_int(end+1) = all(all(abs(vals - vals_true) < tol));

% all radius 0: segments scalar, greater than 1
I_r0 = interval([1;-2],[1;-2]);
segments = 3;
vals = gridPoints(I_r0,segments);
vals_true = [1; -2];
res_int(end+1) = all(all(abs(vals - vals_true) < tol));

% all radius 0: transposed, segments scalar, greater than 1
segments = 3;
vals = gridPoints(I_r0',segments);
vals_true = [1 -2];
res_int(end+1) = all(all(abs(vals - vals_true) < tol));

% 3. interval matrix
% full-dimensional interval matrix: segments scalar, greater than 1
I = interval([-5 -2; 3 -4],[1 2; 5 2]);
segments = 3;
vals = gridPoints(I,segments);
res_intmat = iscell(vals) && all(size(vals{1}) == size(I));

% full-dimensional interval matrix: segments matrix, greater than 1
I = interval([-5 -2; 3 -4],[1 2; 5 2]);
segments = [3 2; 6 4];
vals = gridPoints(I,segments);
res_intmat(end+1) = iscell(vals) && all(size(vals{1}) == size(I));

% combine results
res = res_empty && all(res_int) && all(res_intmat);

% ------------------------------ END OF CODE ------------------------------
