function han = getfcn(obj,options)
% getfcn - returns the function handle of the continuous function specified
%    by the linear probabilistic system object
%
% Syntax:  
%    han = getfcn(obj,options)
%
% Inputs:
%    obj - linProbSys object
%    options - reachability settings
%
% Outputs:
%    han - function handle
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      06-October-2007 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

function dxdt = f(t,x)
    dxdt = obj.A*x+double(obj.B*options.u);
end

han = @f;

end

%------------- END OF CODE --------------