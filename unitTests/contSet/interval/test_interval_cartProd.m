function res = test_interval_cartProd
% test_interval_cartProd - unit test function of Cartesian product
%
% Syntax:  
%    res = test_interval_cartProd
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: see below
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      19-Sep-2019
% Last update:  24-Sep-2019
% Last revision:---

%------------- BEGIN CODE --------------

% TEST 1: Analytical ------------------------------------------------------
% 1. vertcat:
% create intervals
Int1 = interval([-2; -4; -3], ...
                [2;  6;  1]);
Int2 = interval([-1; -5], ...
                [7;  9]);

% compute cartesian product
Int_cartProd = cartProd(Int1, Int2);

% true product
lower_true = [-2; -4; -3; -1; -5];
upper_true = [ 2;  6;  1;  7;  9];
Int_cartProd_true = interval(lower_true, upper_true);

% compare results
res_vertcat = isequal(Int_cartProd,Int_cartProd_true);

% 2. horzcat:
% create intervals
Int1 = interval([-1, -5], ...
                [7,  9]);
Int2 = interval([-2, -4, -3], ...
                [2,  6,  1]);

% compute cartesian product
Int_cartProd = cartProd(Int1, Int2);

% true product
lower_true = [-1,-5,-2,-4,-3];
upper_true = [ 7, 9, 2, 6, 1];
Int_cartProd_true = interval(lower_true, upper_true);

% compare results
res_horzcat = isequal(Int_cartProd,Int_cartProd_true);

% 3. empty interval
% create intervals
Int1 = interval();
Int2 = interval([-1; -5], ...
                [7;  9]);

% compute cartesian product
Int_cartProd = cartProd(Int1, Int2);

% true product
lower_true = [-1; -5];
upper_true = [ 7;  9];
Int_cartProd_true = interval(lower_true, upper_true);

% compare results
res_empty = isequal(Int_cartProd,Int_cartProd_true);

res_analytical = res_vertcat && res_horzcat && res_empty;
% -------------------------------------------------------------------------

% TEST 2: Random ----------------------------------------------------------
% create random interval
dim1 = floor(1 + 9*rand(1));
lower1 = -10*rand(dim1,1);
upper1 = 10*rand(dim1,1);
Int1 = interval(lower1, upper1);
dim2 = floor(1 + 9*rand(1));
lower2 = -10*rand(dim2,1);
upper2 = 10*rand(dim2,1);
Int2 = interval(lower2, upper2);

% compute cartesian product
Int_cartProd = cartProd(Int1, Int2);

% true product
Int_cartProd_true = interval([lower1; lower2], [upper1; upper2]);

% compare results
res_rand = isequal(Int_cartProd,Int_cartProd_true);
% -------------------------------------------------------------------------

% add results
res = res_analytical && res_rand;

if res
    disp('test_cartProd successful');
else
    disp('test_cartProd failed');
end

%------------- END OF CODE --------------