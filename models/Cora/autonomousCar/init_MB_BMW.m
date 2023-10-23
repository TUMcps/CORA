function x0 = init_MB_BMW(initState)
% init_MB - generates the initial state vector for the multi-body model
%
% Syntax:  
%     x0 = init_MB(initState, p)
%
% Inputs:
%     initState - core initial states
%     p - parameter vector
%
% Outputs:
%     x0 - initial state vector
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      11-January-2017
% Last update:  15-December-2017
%               16-June-2023
% Last revision:---

%------------- BEGIN CODE --------------

%get model paramters
p = BMWparameters();

%states
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

%u1 = v_delta steering angle velocity of front wheels
%u2 = acceleration


%obtain initial states from vector
sx0 = initState(1);
sy0 = initState(2);
delta0 = initState(3);
vel0 = initState(4);
Psi0 = initState(5);
dotPsi0 = initState(6);
beta0 = initState(7);


%create equivalent bicycle parameters
g = 9.81; %[m/s^2]

%auxiliary initial states
F0_z_f = p.m_s*g*p.b/((p.a + p.b)) + p.m_uf*g;
F0_z_r = p.m_s*g*p.a/((p.a + p.b)) + p.m_ur*g;

%sprung mass states
x0(1,1) = sx0; % s_x x-position in a global coordinate system
x0(2,1) = sy0; % s_y y-position in a global coordinate system
x0(3,1) = delta0; % steering angle of front wheels
x0(4,1) = cos(beta0)*vel0; % u velocity in x-direction
x0(5,1) = Psi0; % Ψ yaw angle
x0(6,1) = dotPsi0; % Ψ yaw rate
x0(7,1) = 0; % ΦS roll angle
x0(8,1) = 0; % ΦS roll rate
x0(9,1) = 0; % ΘS pitch angle
x0(10,1) = 0; % ΘS pitch rate
x0(11,1) = sin(beta0)*vel0; % v velocity in y-direction
x0(12,1) = 0; % zS z-position (zero height corresponds to steady state solution)
x0(13,1) = 0; % w velocity in z-direction

%unsprung mass states (front)
x0(14,1) = 0; % ΦUF roll angle front
x0(15,1) = 0; % ΦUF roll rate front
x0(16,1) = sin(beta0)*vel0 + p.a*dotPsi0; % vUF velocity in y-direction front
x0(17,1) = (F0_z_f)/(2*p.K_zt); % zUF z-position front
x0(18,1) = 0; % wUF velocity in z-direction front

%unsprung mass states (rear)
x0(19,1) = 0; % ΦUR roll angle rear
x0(20,1) = 0; % ΦUR roll rate rear
x0(21,1) = sin(beta0)*vel0 - p.b*dotPsi0; % vUR velocity in y-direction rear
x0(22,1) = (F0_z_r)/(2*p.K_zt); % zUR z-position rear
x0(23,1) = 0; % wUR velocity in z-direction rear

%wheel states
x0(24,1) = x0(4)/(p.R_w); % ωLF left front wheel angular speed
x0(25,1) = x0(4)/(p.R_w); % ωRF right front wheel angular speed
x0(26,1) = x0(4)/(p.R_w); % ωLR left rear wheel angular speed
x0(27,1) = x0(4)/(p.R_w); % ωRR right rear wheel angular speed

x0(28,1) = 0; % delta_y_f
x0(29,1) = 0; % delta_y_r



%------------- END OF CODE --------------