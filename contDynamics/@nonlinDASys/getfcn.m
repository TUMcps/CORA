function han = getfcn(obj,options)
% getfcn - returns the function handle of the continuous function specified
%    by the DAE system object
%
% Syntax:
%    han = getfcn(obj,options)
%
% Inputs:
%    obj - nonlinDASys object
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
% Written:       17-November-2011
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    function dxdt = f(t, z)
        %obtain x and y
        x = z(1:obj.dim);
        y = z((obj.dim+1):(obj.dim+obj.nrOfConstraints));
        
        %return derivatives
        dxdt(1:obj.dim,1) = obj.dynFile(x, y, options.u);
        dxdt((obj.dim+1):(obj.dim+obj.nrOfConstraints),1) = obj.conFile(x, y, options.u);
    end

    han = @f;

end

% ------------------------------ END OF CODE ------------------------------
