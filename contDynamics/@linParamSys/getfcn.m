function [handle] = getfcn(obj,options)
% getfcn - returns the function handle of the continuous function specified
% by the linear system object
%
% Syntax:  
%    [handle] = getfcn(obj)
%
% Inputs:
%    obj - linIntSys object
%
% Outputs:
%    handle - function handle
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 07-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

function dxdt = f(t,x)
    dxdt = obj.sampleMatrix.A*x+double(obj.B*options.u);
end

handle = @f;
end

%------------- END OF CODE --------------