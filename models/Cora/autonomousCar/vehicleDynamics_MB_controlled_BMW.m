function f = vehicleDynamics_MB_controlled_BMW(x,u)
% vehicleDynamics_MB_controlled_BMW - enhances multi-body model with 
% control for trajectory tracking
%
% Syntax:  
%    f = vehicleDynamics_MB_controlled_BMW(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector (here: reference trajectory)
%
% Outputs:
%    f - time-derivative of the state vector
%
% References:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: DOTcontrol, vehicleDynamics_MB_BMW

% Author:       Matthias Althoff
% Written:      26-August-2011
% Last update:  15-June-2023
% Last revision:---

%------------- BEGIN CODE --------------

%obtain control inputs
carInput = DOTcontrol(x,u);

%simulate vehicle dynamics
f = vehicleDynamics_MB_BMW(x,carInput);



%------------- END OF CODE --------------