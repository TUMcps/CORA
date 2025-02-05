function han = getfcn(nlnsysDA,params)
% getfcn - returns the function handle of the continuous function specified
%    by the DAE system object
%
% Syntax:
%    han = getfcn(nlnsysDA,params)
%
% Inputs:
%    nlnsysDA - nonlinDASys object
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
% Written:       17-November-2011
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    function dxdt = f(t, z)
        %obtain x and y
        x = z(1:nlnsysDA.nrOfDims);
        y = z((nlnsysDA.nrOfDims+1):(nlnsysDA.nrOfDims+nlnsysDA.nrOfConstraints));
        
        %return derivatives
        dxdt(1:nlnsysDA.nrOfDims,1) = nlnsysDA.dynFile(x, y, params.u);
        dxdt((nlnsysDA.nrOfDims+1):(nlnsysDA.nrOfDims+nlnsysDA.nrOfConstraints),1) = nlnsysDA.conFile(x, y, params.u);
    end
    
    han = @f;

end

% ------------------------------ END OF CODE ------------------------------
