function f = DOTBicycleDynamics_controlled_SRX_velEq(t,x,u)
% DOTBicycleDynamics_controlled_SRX_vel - enhances bicycle model with control for
% trajectory tracking
%
% Syntax:  
%    f = DOTBicycleDynamics_controlled_SRX_vel(t,x,u)
%
% Inputs:
%    t - time
%    x - state
%    u - reference trajectory
%    y - sensor noise
%
% Outputs:
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      01-March-2012
% Last update:  15-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%obtain control inputs
carInput = DOTcontrol_SRX_velEq(x,u);

%simulate vehicle dynamics
f = DOTBicycleDynamics_SRX_velEq(t,x,carInput);



%------------- END OF CODE --------------
