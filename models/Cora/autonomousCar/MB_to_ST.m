function x0 = MB_to_ST(initState)
% MB_to_ST - converts the initial state vector of the multi-body model
% to the single-track model
%
% Syntax:  
%     x0 = MB_to_ST(initState)
%
% Inputs:
%     initState - initial states of multi-body model
%
% Outputs:
%     x0 - initial state vector of single-track model
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      16-June-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%states of multi-body model
%x1 = s_x x-position in a global coordinate system
%x2 = s_y y-position in a global coordinate system
%x3 = δ steering angle of front wheels
%x4 = u velocity in x-direction
%x5 = Ψ yaw angle
%x6 = Ψ yaw rate

%x7 = ΦS roll angle
%x8 = ΦS roll rate
%x9 = ΘS pitch angle
%x10 = ΘS pitch rate
%x11 = v velocity in y-direction
%x12 = zS z-position
%x13 = w velocity in z-direction

%x14 = ΦUF roll angle front
%x15 = ΦUF roll rate front
%x16 = vUF velocity in y-direction front
%x17 = zUF z-position front
%x18 = wUF velocity in z-direction front

%x19 = ΦUR roll angle rear
%x20 = ΦUR roll rate rear
%x21 = vUR velocity in y-direction rear
%x22 = zUR z-position rear
%x23 = wUR velocity in z-direction rear

%x24 = ωLF left front wheel angular speed
%x25 = ωRF right front wheel angular speed
%x26 = ωLR left rear wheel angular speed
%x27 = ωRR right rear wheel angular speed

%x28 = delta_y_f
%x29 = delta_y_r

%states of single-track model
%x1 = s_x x-position in a global coordinate system
%x2 = s_y y-position in a global coordinate system
%x3 = δ steering angle of front wheels
%x4 = u velocity in x-direction
%x5 = Ψ yaw angle
%x6 = Ψ yaw rate
%x7 = β slip angle at vehicle center


x0(1,1) = initState(1); % x-position
x0(2,1) = initState(2); % y-position
x0(3,1) = initState(3); % steering angle
x0(4,1) = sqrt(initState(4).^2 + initState(11).^2); % velocity in x-direction
x0(5,1) = initState(5); % yaw angle
x0(6,1) = initState(6); % yaw rate
x0(7,1) = atan(initState(11)./initState(4)); % slip angle

%------------- END OF CODE --------------