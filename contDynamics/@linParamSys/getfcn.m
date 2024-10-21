function han = getfcn(sys,params)
% getfcn - returns the function handle of the continuous function specified
%    by the linear system object
%
% Syntax:
%    han = getfcn(sys)
%
% Inputs:
%    sys - linParamSys object
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
% Written:       07-May-2007 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

han = @(t,x) sys.sampleMatrix.A*x+double(sys.B*params.u);

% ------------------------------ END OF CODE ------------------------------
