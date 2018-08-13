function obj = intersectionTimeParameters(obj,A,u,x_c,x_0_set)
% intersectionTimeParameters - finds the parameters of a polynomial model 
% to determine the time when a linear system hits an arbitrary switching 
% surface h(x)
%
% Syntax:  
%    intersectionTimeParameters(obj,A,u,x_0_set)
%
% Inputs:
%    obj - switchingSurface object
%    A - system matrix
%    u - constant input
%    x_0_set - set of initial points
%
% Outputs:
%
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      21-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%init
y0 = zeros(obj.nrOfParameters,1); % same as linear solution;
converged = 0;

while ~converged
    
    for iParam = 1:obj.nrOfParameters
        h(iParam,1) = obj.timeModelFile(A,u,x_c,x_0_set{iParam},y0);
        J(iParam,:) = obj.jacobianFile(A,u,x_c,x_0_set{iParam},y0);
    end

    %evaluate jacobian
    delta_y = -J\h;
    
    %check convergence
    if norm(delta_y)<1e-10
        converged = 1;
    end
    
    %update steady state solution
    y0 = y0 + delta_y;
end
%set y value
obj.y = y0;



%------------- END OF CODE --------------