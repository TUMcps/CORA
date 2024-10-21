function res = test_isequalFunctionHandle
% test_isequalFunctionHandle - unit test function for equality check of
%    two function handles
%
% Syntax:
%    res = test_isequalFunctionHandle()
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
% Written:       07-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check against itself
% (note: this also holds true for isequal(...))
f = @(x) x;
assert(isequalFunctionHandle(f,f));

% same function, different provinence
f = @(x,u) [x(1) - u(1); x(1)*x(2)];
g = @(x,u) [x(1) - u(1); x(1)*x(2)];
assert(isequalFunctionHandle(f,g));

% different ordering of expressions
f = @(x,u) [x(1) - u(1); x(1)*x(2)];
g = @(x,u) [-u(1) + x(1); x(2)*x(1)];
assert(isequalFunctionHandle(f,g));

% slightly different functions
f = @(x,u) [x(1) - u(1); x(1)*x(2)];
g = @(x,u) [x(1) - u(1); x(1)/x(2)];
assert(~isequalFunctionHandle(f,g));

% different number of input arguments
f = @(x,u) x(1) - u(1);
g = @(x) x(1);
assert(~isequalFunctionHandle(f,g));

% different number of output arguments
f = @(x,u) x(1) - u(1);
g = @(x,u) [x(1) - u(1); x(1)];
assert(~isequalFunctionHandle(f,g));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
