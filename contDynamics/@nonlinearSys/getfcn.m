function han = getfcn(obj,options)
% getfcn - returns the function handle of the continuous function specified
%    by the linear system object
%
% Syntax:
%    han = getfcn(obj)
%
% Inputs:
%    obj - nonlinearSys object
%    options - reachability settings
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
% Written:       17-October-2007 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    function dxdt = f(t,x)
        dxdt = obj.mFile(x,options.u);
    end
    
    han = @f;

end

% ------------------------------ END OF CODE ------------------------------
