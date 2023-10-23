function f = vehicleDynamics_ST_controlled_BMW(x,u)
% vehicleDynamics_ST_controlled_BMW - enhances single-track model with 
% control for trajectory tracking; adapted from [1] and changed according
% to vehicle model in CommonRoad
%
% Syntax:  
%    f = vehicleDynamics_ST_controlled_BMW(x,u)
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
%               21-June-2023
% Last revision:---

%------------- BEGIN CODE --------------

% obtain control inputs
carInput = DOTcontrol(x,u);

% obtain disturbances
w = u(11:end);

% combine control input and disturbance
u = [carInput; w];

% simulate vehicle dynamics
f = vehicleDynamics_ST_BMW(x,u);



%------------- END OF CODE --------------