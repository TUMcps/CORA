function HA = gearbox()
% gearbox - gearbox benchmark from [1], which is part of ARCH competition
%
% Syntax:  
%    HA = gearbox()
%
% Inputs:
%    -
%
% Outputs:
%    HA - hybrid automaton
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: 
%   [1] H. Chen and et al. "Motor-transmission drive system: a benchmark 
%       example for safety verification"
%
% See also: example_hybrid_reach_ARCH20_gearbox_GRBX01

% Author:        Matthias Althoff, Niklas Kochdumper
% Written:       19-Dec-2016 
% Last update:   13-May-2020 (NK, moved to seperate file)
% Last revision: ---

%------------- BEGIN CODE --------------


% Parameter ---------------------------------------------------------------

theta = 0.9;                    % coefficient of restitution
m_s = 3.2;                      % [kg], mass of sleve
m_g2 = 18.1;                    % [kg], mass of second gear
J_g2 = 0.7;                     % [kg m^2], inertia of second gear
R_s = 0.08;                     % [m], radius of sleeve
Theta = 36/180*pi;              % [rad], included angle of gear
b = 0.01;                       % [m], width of gear spline
delta_p = -0.003;               % [m], p_x sleeve meshes with gear 
n = 0;                          % integer number in guard
F_s = 70;                       % [N], shifting force
T_f = 1;                        % [Mm], resisting moment




% Continuous Dynamics -----------------------------------------------------

% system matrix
A = [0 0 0 0 0;...
     0 0 0 0 0;...
     1 0 0 0 0;...
     0 1 0 0 0;...
     0 0 0 0 0]; 
    
B = [1/m_s; 0; 0; 0; 0];

c = B*F_s + [0; -R_s*T_f/J_g2; 0; 0; 0];

% get dimension
dim_x = length(A);

% instantiate linear dynamics 
linSys  = linearSys('linearSys',A, B, c);




% Transition 1 ------------------------------------------------------------

%define large and small distance
dist = 1e3;
smallDist = 1e3*eps;

% guard set
ch_n = [0, 0, -tan(Theta), -1, 0];
ch_d = -2*n*b;
ch_C(1,:) = [0, 0, 1, 0, 0];
ch_C(2,:) = [0, 0, -1, 0, 0];
ch_C(3,:) = [-sin(Theta), -cos(Theta), 0, 0, 0];
ch_D = [delta_p; b/tan(Theta); 0];
guard1 = conHyperplane(halfspace(ch_n,ch_d),ch_C,ch_D);

% inequality constraints for invariant set
C1(1,:) = [0, 0, -tan(Theta), -1, 0];
C1(2,:) = [0, 0, 1, 0, 0];
C1(3,:) = [0, 0, -1, 0, 0];
C1(4,:) = [-sin(Theta), -cos(Theta), 0, 0, 0];
C1(5,:) = -C1(1,:);
d1 = [-2*n*b; delta_p; b/tan(Theta); 0; -2*n*b + smallDist];

% reset function
denominator = m_s*cos(Theta)^2 + m_g2*sin(Theta)^2;
a_11 = (m_s*cos(Theta)^2 - m_g2*theta*sin(Theta)^2)/denominator;
a_12 = (-1*(theta+1)*m_g2*sin(Theta)*cos(Theta))/denominator;
a_21 = (-1*(theta+1)*m_s*sin(Theta)*cos(Theta))/denominator;
a_22 = (m_g2*sin(Theta)^2 - m_s*theta*cos(Theta)^2)/denominator;
a_51 = ((theta+1)*m_s*m_g2*sin(Theta))/denominator;
a_52 = ((theta+1)*m_s*m_g2*cos(Theta))/denominator;
reset1.A = [a_11, a_12, 0, 0, 0; ...
    a_21, a_22, 0, 0, 0; ...]
    0, 0, 1, 0, 0; ...
    0, 0, 0, 1, 0; ...
    a_51, a_52, 0, 0, 1]; 
reset1.c = zeros(dim_x,1);

% transition 
trans(1) = transition(guard1,reset1,1);



% Transition 2 ------------------------------------------------------------

% guard set
ch_n = [0, 0, -tan(Theta), 1, 0];
ch_d = 2*n*b;
ch_C(1,:) = [0, 0, 1, 0, 0];
ch_C(2,:) = [0, 0, -1, 0, 0];
ch_C(3,:) = [-sin(Theta), +cos(Theta), 0, 0, 0];
ch_D = [delta_p; b/tan(Theta); 0];
guard2 = conHyperplane(halfspace(ch_n,ch_d),ch_C,ch_D);

% inequality constraints for invariant set
C2(1,:) = [0, 0, -tan(Theta), 1, 0];
C2(2,:) = [0, 0, 1, 0, 0];
C2(3,:) = [0, 0, -1, 0, 0];
C2(4,:) = [-sin(Theta), +cos(Theta), 0, 0, 0];
C2(5,:) = -C2(1,:);
d2 = [2*n*b; delta_p; b/tan(Theta); 0; 2*n*b + smallDist];

% reset function
denominator = m_s*cos(Theta)^2 + m_g2*sin(Theta)^2;
a_11 = (m_s*cos(Theta)^2 - m_g2*theta*sin(Theta)^2)/denominator;
a_12 = ((theta+1)*m_g2*sin(Theta)*cos(Theta))/denominator;
a_21 = ((theta+1)*m_s*sin(Theta)*cos(Theta))/denominator;
a_22 = (m_g2*sin(Theta)^2 - m_s*theta*cos(Theta)^2)/denominator;
a_51 = ((theta+1)*m_s*m_g2*sin(Theta))/denominator;
a_52 = (-1*(theta+1)*m_s*m_g2*cos(Theta))/denominator;
reset2.A = [a_11, a_12, 0, 0, 0; ...
    a_21, a_22, 0, 0, 0; ...]
    0, 0, 1, 0, 0; ...
    0, 0, 0, 1, 0; ...
    a_51, a_52, 0, 0, 1]; 
reset2.c = zeros(dim_x,1);

% transition
trans(2) = transition(guard2,reset2,1);



% Transition 3 ------------------------------------------------------------

% guard set
C3(1,:) = [0, 0, -1, 0, 0];
d3 = [-delta_p]; % different from paper!!
guard3 = conHyperplane(C3,d3);

% reset function
reset3.A = [0, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0; ...
    0, 0, 1, 0, 0; ...
    0, 0, 0, 1, 0; ...
    m_s, m_s, 0, 0, 1]; 
reset3.c = zeros(dim_x,1);

% transition
trans(3) = transition(guard3,reset3,2);



% Hybrid Automaton --------------------------------------------------------

% invariant set
Cinv = [-C1(1,:); -C2(1,:); -C1(3,:)];
dinv = [d1(1); d2(1); d1(3)];
dinv = dinv + ones(length(dinv),1)*1e-3;
inv = polytope(Cinv,dinv);

% location
loc(1) = location('loc1',inv,trans,linSys);

loc(2) = location('loc2',inv,transition(),linearSys(zeros(5),zeros(5,1)));

% hybrid automaton
HA = hybridAutomaton(loc);

%------------- END OF CODE --------------