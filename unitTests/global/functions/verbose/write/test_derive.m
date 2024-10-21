function res = test_derive
% test_derive - unit tests for the derivation of symbolic functions
%
% Syntax:
%    res = test_derive
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
fname = 'test_derive_generatedfile';
% fullname needs to match generated file, otherwise delete won't work
fullname = [path filesep fname '.m'];
% to delete generated file: delete(fullname); 

% 1D -> 1D, univariate
x = sym('x',[1,1]);
f = @(x) x^3;
f_sym = x(1)^3;
J_han = derive('FunctionHandle',f,'Vars',{x},'Path','none');
J_sym = derive('SymbolicFunction',f_sym,'Vars',{x},'Path','none');
J_true = 3*x(1)^2;
assert(numel(J_han) == 1);
assert(isequal(J_han{1},J_true));
assert(numel(J_sym) == 1);
assert(isequal(J_sym{1},J_true));

% 1D -> 1D, multivariate
x = sym('x',[1,1]);
u = sym('u',[1,1]);
f = @(x,u) x(1)^3*u(1);
f_sym = x(1)^3*u(1);
J_han = derive('FunctionHandle',f,'Vars',{x,u},'Path','none');
J_sym = derive('SymbolicFunction',f_sym,'Vars',{x,u},'Path','none');
J_true = {3*u(1)*x(1)^2, x(1)^3};
assert(numel(J_han) == 2);
assert(isequal(J_han{1},J_true{1}));
assert(isequal(J_han{2},J_true{2}));
assert(numel(J_sym) == 2);
assert(isequal(J_sym{1},J_true{1}));
assert(isequal(J_sym{2},J_true{2}));

% 2D -> 2D, univariate
x = sym('x',[2,1]);
u = sym('u',[1,1]);
f = @(x) [x(1)*x(2), -2*x(1)^3];
f_sym = [x(1)*x(2); -2*x(1)^3];
J_han = derive('FunctionHandle',f,'Vars',{x},'Path','none');
J_sym = derive('SymbolicFunction',f_sym,'Vars',{x},'Path','none');
J_true = {[x(2), x(1); -6*x(1)^2, 0]};
assert(numel(J_han) == 1);
assert(isequal(J_han{1},J_true{1}));
assert(numel(J_sym) == 1);
assert(isequal(J_sym{1},J_true{1}));

% 2D -> 2D, multivariate
x = sym('x',[2,1]);
u = sym('u',[1,1]);
f = @(x,u) [x(1)*x(2) - u(1); -2*x(1)^3 + x(2)*u(1)];
f_sym = [x(1)*x(2) - u(1); -2*x(1)^3 + x(2)*u(1)];
J_han = derive('FunctionHandle',f,'Vars',{x,u},'Path','none');
J_sym = derive('SymbolicFunction',f_sym,'Vars',{x,u},'Path','none');
J_true = {[x(2), x(1); -6*x(1)^2, u(1)], [-1; x(2)]};
assert(numel(J_han) == 2);
assert(isequal(J_han{1},J_true{1}));
assert(isequal(J_han{2},J_true{2}));
assert(numel(J_sym) == 2);
assert(isequal(J_sym{1},J_true{1}));
assert(isequal(J_sym{2},J_true{2}));

% 2Dx2D -> 2Dx2D, univariate
x = sym('x',[2,1]);
f = @(x) [x(1)*x(2), -x(1)^3; x(2)^2-x(1), -x(2)*x(1)];
f_sym = [x(1)*x(2), -x(1)^3; x(2)^2-x(1), -x(2)*x(1)];
J_han = derive('FunctionHandle',f,'Vars',{x},'Path','none');
J_sym = derive('SymbolicFunction',f_sym,'Vars',{x},'Path','none');
J_true = {sym(zeros(2,2,2))};
J_true{1}(:,:,1) = [x(2), -3*x(1)^2; -1, -x(2)];
J_true{1}(:,:,2) = [x(1), 0; 2*x(2), -x(1)];
assert(numel(J_han) == 1);
assert(isequal(J_han{1},J_true{1}));
assert(numel(J_sym) == 1);
assert(isequal(J_sym{1},J_true{1}));

% 2Dx2D -> 2Dx2D, multivariate
x = sym('x',[2,1]);
u = sym('u',[1,1]);
f = @(x,u) [x(1)*x(2) - u(1), -x(1)^3; x(2)^2*u(1)-x(1), -x(2)*x(1)];
f_sym = [x(1)*x(2) - u(1), -x(1)^3; x(2)^2*u(1)-x(1), -x(2)*x(1)];
J_han = derive('FunctionHandle',f,'Vars',{x,u},'Path','none');
J_sym = derive('SymbolicFunction',f_sym,'Vars',{x,u},'Path','none');
J_true = {sym(zeros(2,2,2))};
J_true{1}(:,:,1) = [x(2), -3*x(1)^2; -1, -x(2)];
J_true{1}(:,:,2) = [x(1), 0; 2*x(2)*u(1), -x(1)];
J_true{2}(:,:,1) = [-1, 0; x(2)^2, 0];
assert(numel(J_han) == 2);
assert(isequal(J_han{1},J_true{1}));
assert(isequal(J_han{2},J_true{2}));
assert(numel(J_sym) == 2);
assert(isequal(J_sym{1},J_true{1}));
assert(isequal(J_sym{2},J_true{2}));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
