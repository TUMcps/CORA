function [p,omega_ref_center,omega_ref_delta,x1_0,delta_x1_0,xy_0,delta_xy_0,delta_T_m0] = powertrainParameters()

%set parameters of the powertrain
p.alpha = 0.03;     %backlash size (half gap) [rad]
p.tau_eng = 0.1;    %time constant of the engine [s]
p.b_l = 5.6;        %viscous friction of wheels [Nm/(rad/s)]
p.b_m = 0;          %viscous friction of engine [Nm/(rad/s)]
p.b_i = 1;          %viscous friction of additional inertias [Nm/(rad/s)]
p.i = 12;           %transmission ratio, Theta_m/Theta_1 [rad/rad]
p.k = 10e3;         %shaft stiffness [Nm/rad]
p.k_i = 10e4;       %shaft stiffness of additional inertias [Nm/rad]
%p.k_i = 10e3;      %shaft stiffness of additional inertias [Nm/rad]
p.r = 0.33;         %wheel radius [m]
p.J_l = 140;        %moment of inertia of wheels and vehicle mass [kgm^2]
p.J_m = 0.3;        %moment of inertia of engine flywheel [kgm^2]
p.J_i = 0.01;       %moment of inertia of additional inertias [kgm^2]

%control parameters:
p.k_P = 0;
p.k_I = 0;
p.k_D = 0;
p.k_K = 0.5;
p.k_KD = 0.5;
p.k_KI = 0.5;
p.k_KK = 0;

%create initial zonotope
omega_ref_center = 30;
omega_ref_delta = 10;
T_l_center = -300;

%initial set for loc 3:
x1_0 = p.b_l/p.k*omega_ref_center - p.alpha + T_l_center/p.k;
delta_x1_0 = p.b_l/p.k*omega_ref_delta;
xy_0 = p.b_l/p.k_i*omega_ref_center + T_l_center/p.k_i;
delta_xy_0 = p.b_l/p.k_i*omega_ref_delta;


p.T_m0 = p.k/p.i*(x1_0 + p.alpha) + p.b_m*omega_ref_center; %needs to be integrated in U
delta_T_m0 = p.k/p.i*delta_x1_0 + p.b_m*omega_ref_delta;    %needs to be integrated in U

end
