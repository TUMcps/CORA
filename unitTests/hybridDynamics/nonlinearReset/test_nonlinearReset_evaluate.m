function res = test_nonlinearReset_evaluate
% test_nonlinearReset_evaluate - test function for evaluation of nonlinear
%    reset functions
%
% Syntax:
%    res = test_nonlinearReset_evaluate
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
% Written:       13-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-14;

% set path to this folder
path = mfilename('fullpath');
idxFilesep = strfind(path,filesep);
path = path(1:idxFilesep(end));
fname = 'test_nonlinearReset_evaluate_generatedfile';
% fullname needs to match generated file, otherwise delete won't work
fullname_jacobian = [path filesep fname '_jacobian.m'];
fullname_hessian = [path filesep fname '_hessian.m'];
% fullname_thirdorder = [path filesep fname 'thirdOrderTensor_.m'];

% only states
f = @(x,u) [x(1)*x(2); x(2)];
nonlinReset = nonlinearReset(f);
nonlinReset = derivatives(nonlinReset,path,fname,2);
x = [1;-2];
x_ = evaluate(nonlinReset,x);
x_true = f(x);
assert(compareMatrices(x_,x_true,tol,"equal",true));

% clean up
delete(fullname_jacobian);
delete(fullname_hessian);


% states and inputs
% f = @(x,u) [x(1)^3*x(2)^4*u(1)^5*u(2)^6; x(1)^4*x(2)^5*u(1)^6*u(2)^3];
f = @(x,u) [x(1)*x(2)^2-u(1); -x(1)^3 + x(2)*u(2)];
nonlinReset = nonlinearReset(f);
nonlinReset = derivatives(nonlinReset,path,fname,2);

% point-wise evaluation
x = [1;-2]; u = [5;4];
x_ = evaluate(nonlinReset,x,u);
x_true = f(x,u);
assert(compareMatrices(x_,x_true,tol,"equal",true));

% set-based evaluation
x = zonotope([100;-2],[0.1 0.2; -0.4 0]);
u = zonotope([5;4],[0.3 -0.1; 0.2 0.1]);
x_ = evaluate(nonlinReset,x,u);
% ...sample random points within x and u and check for containment
numPoints = 100;
x_rand = randPoint(x,numPoints);
u_rand = randPoint(u,numPoints);
x_sampled = zeros(2,numPoints);
for i=1:numPoints
    x_sampled(:,i) = f(x_rand(:,i),u_rand(:,i));
end
assert(all(contains_(x_,x_sampled,'exact',tol)));

% clean up
delete(fullname_jacobian);
delete(fullname_hessian);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
