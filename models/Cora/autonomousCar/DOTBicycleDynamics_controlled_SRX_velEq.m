function f = DOTBicycleDynamics_controlled_SRX_velEq(x,u)
% DOTBicycleDynamics_controlled_SRX_vel - enhances bicycle model (see [1])
%                                         with control for trajectory 
%                                         tracking
%
% Syntax:  
%    f = DOTBicycleDynamics_controlled_SRX_vel(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector (here: reference trajectory)
%
% Outputs:
%    f - time-derivative of the state vector
%
% References:
%    [1] M. Althoff and J. M. Dolan. Online verification of automated
%        road vehicles using reachability analysis.
%        IEEE Transactions on Robotics, 30(4):903-918, 2014.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: DOTcontrol_SRX_velEq, DOTBicycleDynamics_SRX_velEq

% Author:       Matthias Althoff
% Written:      01-March-2012
% Last update:  15-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

    % obtain control inputs
    carInput = DOTcontrol_SRX_velEq(x,u);

    % simulate vehicle dynamics
    f = DOTBicycleDynamics_SRX_velEq(x,carInput);

%------------- END OF CODE --------------