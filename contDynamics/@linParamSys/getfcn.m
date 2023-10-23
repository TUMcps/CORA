function han = getfcn(obj,options)
% getfcn - returns the function handle of the continuous function specified
%    by the linear system object
%
% Syntax:
%    han = getfcn(obj)
%
% Inputs:
%    obj - linParamSys object
%
% Outputs:
%    han - function handle
%
% Example: 
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

han = @(t,x) obj.sampleMatrix.A*x+double(obj.B*options.u);

end

% ------------------------------ END OF CODE ------------------------------
