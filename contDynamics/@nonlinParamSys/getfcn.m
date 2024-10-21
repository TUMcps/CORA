function han = getfcn(sys,params)
% getfcn - returns the function handle of the continuous function specified
%    by the nonlinear parametric system
%
% Syntax:
%    han = getfcn(sys,params)
%
% Inputs:
%    sys - linearSys object
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
% Written:       17-October-2007 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    function dxdt = f(t,x)
        dxdt = sys.mFile(x,params.u,params.p);
    end
    
    han = @f;

end

% ------------------------------ END OF CODE ------------------------------
