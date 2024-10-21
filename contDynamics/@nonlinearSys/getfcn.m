function han = getfcn(nlnsys,params)
% getfcn - returns the function handle of the continuous function specified
%    by the linear system object
%
% Syntax:
%    han = getfcn(nlnsys)
%
% Inputs:
%    nlnsys - nonlinearSys object
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
        dxdt = nlnsys.mFile(x,params.u);
    end
    
    han = @f;

end

% ------------------------------ END OF CODE ------------------------------
