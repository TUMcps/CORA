function nonlinReset = nonlinearReset(linReset)
% nonlinearReset - converts a linearReset object into a nonlinearReset
%    object
%
% Syntax:
%    nonlinReset = nonlinearReset(linReset)
%
% Inputs:
%    linReset - linearReset object
%
% Outputs:
%    nonlinReset - nonlinearReset object
%
% Example:
%    A = [1 0; 0 1]; B = [1; 0]; c = [-1; 1];
%    linReset = linearReset(A,B,c);
%    nonlinReset = nonlinearReset(linReset);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reset, nonlinearReset

% Authors:       Mark Wetzlinger
% Written:       08-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
if linReset.preStateDim == 0
    nonlinReset = nonlinearReset(@(x,u) []);
    return
end

% dirty hack... since we require that the input arguments to the generated
% MATLAB function use in(:,1) for indexing instead of just x1 (which
% results in 'any dimension' as far as inputArgsLength is concerned...)
x_sym = sym('x',[linReset.preStateDim+1,1],'real');
u_sym = sym('u',[linReset.inputDim+1,1],'real');

% evaluate linear reset symbolically, write to anonymous function
f_eval = linReset.A*x_sym(1:end-1) + linReset.B*u_sym(1:end-1) + linReset.c;
if isempty(f_eval)
    f = @(x,u) [];
else
    f = matlabFunction(f_eval,'Vars',{x_sym,u_sym});
end

% instantiate nonlinearReset object
nonlinReset = nonlinearReset(f);

% ------------------------------ END OF CODE ------------------------------
