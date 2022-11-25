function [handle] = getfcn(obj,options)
% getfcn - returns the function handle of the continuous function specified
% by the linear system object
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

% Author:       Matthias Althoff
% Written:      07-May-2007 
% Last update:  20-March-2008
%               05-December-2017
% Last revision:---

%------------- BEGIN CODE --------------

if ~isempty(obj.c)
    c = obj.c;
else
    c = 0;
end

function dxdt = f(~,x)
    dxdt = obj.A*x + obj.B*options.u + c;
end

handle = @f;
end

%------------- END OF CODE --------------