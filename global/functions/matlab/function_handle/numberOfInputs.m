function [count,out] = numberOfInputs(f,varargin)
% numberOfInputs - computes the number of inputs of a function handle;
%    this function is kept for legacy reasons -> use 'inputArgsLength'
%    instead (to avoid confusion regarding the function name)
%
% Syntax:
%    [count,out] = numberOfInputs(f)
%    [count,out] = numberOfInputs(f,inpArgs)
%
% Inputs:
%    f - function handle 
%    inpArgs - number of input arguments for the function (max. 26)
%
% Outputs:
%    count - vector storing the length of each input argument
%    out - output dimension of the function handle
%
% Example:
%    f = @(x,u) [x(1)*x(5)^2; sin(x(3)) + u(2)];
%    numberOfInputs(f)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys

% Authors:       Mark Wetzlinger
% Written:       22-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

warning("Function has been renamed to 'inputArgsLength' to avoid confusion.");
[count,out] = inputArgsLength(f,varargin{:});

% ------------------------------ END OF CODE ------------------------------
