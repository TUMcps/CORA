function res = test_nonlinearReset_derivatives
% test_nonlinearReset_derivatives - test function for derivative
%    computation of nonlinear reset functions
%
% Syntax:
%    res = test_nonlinearReset_derivatives
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
% Written:       12-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path to this folder
path = mfilename('fullpath');
idxFilesep = strfind(path,filesep);
path = path(1:idxFilesep(end));
fname = 'test_nonlinearReset_derivatives_generatedfile';
% fullname needs to match generated file, otherwise delete won't work
fullname_jacobian = [path filesep fname '_jacobian.m'];
fullname_hessian = [path filesep fname '_hessian.m'];
fullname_thirdorder = [path filesep fname 'thirdOrderTensor_.m'];

% empty
% nonlinReset = nonlinearReset();
% assert(nonlinReset.preStateDim == 0);
% assert(nonlinReset.inputDim == 1);
% assert(nonlinReset.postStateDim == 0);

% only states
f = @(x,u) [x(1)*x(2); x(2)];
% only Jacobian and Hessian
nonlinReset = nonlinearReset(f);
nonlinReset = derivatives(nonlinReset,path,fname,2);
assert(~isempty(nonlinReset.J));
assert(~isempty(nonlinReset.H));
delete(fullname_jacobian);
delete(fullname_hessian);

% states and inputs
f = @(x,u) [x(1) - u(1); x(1)*x(2)];
nonlinReset = nonlinearReset(f);
nonlinReset = derivatives(nonlinReset,path,fname,2);
assert(~isempty(nonlinReset.J));
assert(~isempty(nonlinReset.H));
delete(fullname_jacobian);
delete(fullname_hessian);

% states and inputs, different output dimension
f = @(x,u) x(1)*x(2) - u(1);
nonlinReset = nonlinearReset(f);
nonlinReset = derivatives(nonlinReset,path,fname,2);
assert(~isempty(nonlinReset.J));
assert(~isempty(nonlinReset.H));
delete(fullname_jacobian);
delete(fullname_hessian);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
