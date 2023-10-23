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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. vertcat:
% create intervals
I1 = interval([-2; -4; -3], ...
              [2;  6;  1]);
I2 = interval([-1; -5], ...
              [7;  9]);

% compute cartesian product
I_ = cartProd(I1, I2);

% true product
lb_true = [-2; -4; -3; -1; -5];
ub_true = [ 2;  6;  1;  7;  9];
I_true = interval(lb_true, ub_true);

% compare results
res_vertcat = isequal(I_,I_true);

% 2. horzcat:
% create intervals
I1 = interval([-1, -5], ...
                [7,  9]);
I2 = interval([-2, -4, -3], ...
                [2,  6,  1]);

% compute cartesian product
I_ = cartProd(I1, I2);

% true product
lb_true = [-1,-5,-2,-4,-3];
ub_true = [ 7, 9, 2, 6, 1];
I_true = interval(lb_true, ub_true);

% compare results
res_horzcat = isequal(I_,I_true);

% 3. interval-numeric case
I1 = interval([-2; -4; -3], ...
              [2;  6;  1]);
num = [2;1];

% compute cartesian product
I_ = cartProd(I1,num);

% true product
lb_true = [-2; -4; -3; 2; 1];
ub_true = [ 2;  6;  1; 2; 1];
I_true = interval(lb_true,ub_true);

% compare results
res_num1 = isequal(I_,I_true);

% compute cartesian product
I_ = cartProd(num,I1);

% true product
lb_true = [2; 1; -2; -4; -3];
ub_true = [2; 1; 2;  6;  1];
I_true = interval(lb_true,ub_true);

% compare results
res_num2 = isequal(I_,I_true);

% % 4. empty interval
% % create intervals
% Int1 = interval();
% Int2 = interval([-1; -5], ...
%                 [7;  9]);
% 
% % compute Cartesian product
% try
%     cartProd(Int1, Int2);
%     res_empty = false;
% catch ME
%     if ~strcmp(ME.identifier,'CORA:notSupported')
%         rethrow(ME);
%     else
%         res_empty = true;
%     end
% end

res = res_vertcat && res_horzcat && res_num1 && res_num2; % && res_empty;

% ------------------------------ END OF CODE ------------------------------
