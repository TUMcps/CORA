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
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-September-2019
% Last update:   24-September-2019
%                03-January-2023 (MW, add interval-numeric cases)
%                03-December-2023 (MW, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1. vertcat
I1 = interval([-2;-4;-3],[2;6;1]);
I2 = interval([-1;-5],[7;9]);
I_cartProd = cartProd(I1,I2);
I_true = interval([-2;-4;-3;-1;-5],[2;6;1;7;9]);
res(end+1,1) = isequal(I_cartProd,I_true);

% 2. horzcat
I1 = interval([-1,-5],[7,9]);
I2 = interval([-2,-4,-3],[2,6,1]);
I_cartProd = cartProd(I1,I2);
I_true = interval([-1,-5,-2,-4,-3],[7,9,2,6,1]);
res(end+1,1) = isequal(I_cartProd,I_true);

% 3. interval-numeric case
I1 = interval([-2;-4;-3],[2;6;1]);
num = [2;1];
% ...I x num
I_cartProd = cartProd(I1,num);
I_true = interval([-2; -4; -3; 2; 1],[2; 6; 1; 2; 1]);
res(end+1,1) = isequal(I_cartProd,I_true);
% ...num x I
I_cartProd = cartProd(num,I1);
I_true = interval([2;1;-2;-4;-3],[2;1;2;6;1]);
res(end+1,1) = isequal(I_cartProd,I_true);

% 4. unbounded case
I1 = interval([-Inf;-2;1],[-2;Inf;Inf]);
I2 = interval([2;-Inf],[Inf;Inf]);
I_cartProd = cartProd(I1,I2);
I_true = interval([-Inf;-2;1;2;-Inf],[-2;Inf;Inf;Inf;Inf]);
res(end+1,1) = isequal(I_cartProd,I_true);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
