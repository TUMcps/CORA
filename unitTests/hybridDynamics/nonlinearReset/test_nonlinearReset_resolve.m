function res = test_nonlinearReset_resolve
% test_nonlinearReset_resolve - test function for the input resolution
%    of a nonlinear reset function
%
% Syntax:
%    res = test_nonlinearReset_resolve
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
% See also: test_linearReset_resolve

% Authors:       Mark Wetzlinger
% Written:       14-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init synchronized nonlinear reset function
f = @(x,u) [x(1) - u(1); ...
            x(2); ...
            x(2) + x(3)*u(2); ...
            x(4) + u(3); ...
            x(5)];
nonlinReset = nonlinearReset(f);

% binds and flows
stateBinds = {[1,2,3],[4,5]};
inputBinds = {[2 1; 0 1],[1 1]};
sys1 = nonlinearSys(@(x,u) [x(1); x(2)+u(1); x(3)-u(2)], ...
                    @(x,u) x(1)*x(2)^2 + sqrt(x(3))*u(2));
sys2 = nonlinearSys(@(x,u) [x(1)+u(1); x(2)], ...
                    @(x,u) -x(1)*x(2));
flowList = {sys1;sys2};

% resolve input binds
nonlinReset_ = resolve(nonlinReset,flowList,stateBinds,inputBinds);
f_res = @(x,u) [x(1) + x(4)*x(5); ...
                x(2); ...
                x(2) + x(3)*u(1); ...
                x(4) + x(1)*x(2)^2 + sqrt(x(3))*u(1); ...
                x(5)];
assert(isequalFunctionHandle(nonlinReset_.f,f_res));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
