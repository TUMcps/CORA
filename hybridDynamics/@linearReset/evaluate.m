function x_ = evaluate(linReset,x,varargin)
% evaluate - evaluates the linear reset function for a given state (set)
%    and input (set)
%
% Syntax:
%    x_ = evaluate(linReset,x)
%    x_ = evaluate(linReset,x,u)
%
% Inputs:
%    nonlinReset - nonlinearReset object
%    x - state before reset
%    u - (optional) input before reset
%
% Outputs:
%    x_ - state after reset
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearReset/evaluate

% Authors:       Mark Wetzlinger
% Written:       07-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default value
narginchk(2,3);
u = setDefaultValues({zeros(linReset.inputDim,1)},varargin);

% evaluate reset function
x_ = linReset.A * x + linReset.B * u + linReset.c;

% ------------------------------ END OF CODE ------------------------------
