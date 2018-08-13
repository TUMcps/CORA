function [handle] = getfcn(obj,options)
% getfcn - returns the function handle of the continuous function specified
% by the zero dynamics system object
%
% Syntax:  
%    [handle] = getfcn(obj)
%
% Inputs:
%    obj - linearSys object
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
% Written: 26-November-2008 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

function dxdt = f(t,x)
    dxdt = 0*x;
end

handle = @f;
end

%------------- END OF CODE --------------