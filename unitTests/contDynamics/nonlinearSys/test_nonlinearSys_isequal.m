function res = test_nonlinearSys_isequal
% test_nonlinearSys_isequal - unit test for equality check
%
% Syntax:
%    res = test_nonlinearSys_isequal
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
% Written:       09-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% instantiate nonlinear functions

% 1D, without inputs
f = @(x,u) x(1)^2;
sys1 = nonlinearSys(f);

% 1D, with inputs
f = @(x,u) x(1)^2 - u(1);
sys2 = nonlinearSys(f);

% 1D, with inputs and output equation
g = @(x,u) x(1)*u(1);
sys3 = nonlinearSys(f,g);

% 2D, without inputs
f = @(x,u) [x(1)^2 - x(2); sqrt(x(2) - x(1))];
sys4 = nonlinearSys(f);

% 2D, with inputs
f = @(x,u) [x(2)^2 - u(1); u(2)*u(3) - x(1)];
sys5 = nonlinearSys(f);
% re-ordered version
f = @(x,u) [-u(1) + x(2)^2; -x(1) + u(2)*u(3)];
sys5_ = nonlinearSys(f);

% 2D, with input and output equation
g = @(x,u) [x(1)*u(1); sqrt(x(2))];
sys6 = nonlinearSys(f,g);


% check results
res = isequal(sys1,sys1);
res(end+1,1) = ~isequal(sys1,sys2);
res(end+1,1) = ~isequal(sys2,sys3);
res(end+1,1) = ~isequal(sys4,sys5);
res(end+1,1) = isequal(sys5,sys5_);
res(end+1,1) = ~isequal(sys5_,sys6);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
