function han = getfcn(sys,params)
% getfcn - returns the function handle of the continuous function specified
%    by the linear probabilistic system object
%
% Syntax:
%    han = getfcn(sys,params)
%
% Inputs:
%    sys - linProbSys object
%    params - model parameters
%
% Outputs:
%    han - function handle
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       06-October-2007 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

han = @(t,x) sys.A*x+double(sys.B*params.u);

% ------------------------------ END OF CODE ------------------------------
