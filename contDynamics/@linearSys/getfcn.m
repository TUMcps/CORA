function han = getfcn(linsys,params)
% getfcn - returns the function handle of the continuous function specified
%    by the linear system object
%
% Syntax:
%    han = getfcn(linsys,params)
%
% Inputs:
%    linsys - linearSys object
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
% Last update:   20-March-2008
%                05-December-2017
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

han = @(t,x) linsys.A*x + linsys.B*params.u + linsys.c + linsys.E*params.w;

% ------------------------------ END OF CODE ------------------------------
