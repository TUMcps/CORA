function f=highorderBicycleDynamics(t,x,u)

% DOTBicycleDynamics - generates bicycle model for comparison with the DOT 
% (department of transportation) car dynamics
%
% Syntax:  
%    f = DOTBicycleDynamics(t,x,u)
%
% Inputs:
%    f
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
% Written:      22-August-2011
% Last update:  30-August-2011
%               18-Januar-2012
%               21-July-2017 M.Klischat: added state for 7th generator
% Last revision:---

%------------- BEGIN CODE --------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%masses
p.m = lb_sec2_ft_IN_kg(74.91452); %vehicle mass [kg]
p.m_s = lb_sec2_ft_IN_kg(66.17221); %sprung mass [kg]
p.m_uf = lb_sec2_ft_IN_kg(4.371153); %unsprung mass front [kg]
p.m_ur = lb_sec2_ft_IN_kg(4.371153); %unsprung mass rear [kg]

%axes distances
p.a = ft_IN_m(3.793293); %distance from spring mass center of gravity to front axle [m]
p.b = ft_IN_m(4.667707); %distance from spring mass center of gravity to rear axle [m]

%moments of inertia of sprung mass
p.I_Phi_s = lb_ft_sec2_IN_kg_m2(152.871); %moment of inertia for sprung mass in roll [kg m^2]
p.I_y_s = lb_ft_sec2_IN_kg_m2(1154.888); %moment of inertia for sprung mass in pitch [kg m^2]
p.I_z = lb_ft_sec2_IN_kg_m2(1321.416); %moment of inertia for sprung mass in yaw [kg m^2]
p.I_xz_s = lb_ft_sec2_IN_kg_m2(0); %moment of inertia cross product [kg m^2]

%steering parameters
p.K_str = 20.8; %overall steering ratio
p.K_scf = rad_ft_lb_IN_rad_sec2_kg_m2(0.00005); %steering compliance for steering gear (front) [rad s^2/(kg m2)] 
p.K_scb = rad_ft_lb_IN_rad_sec2_kg_m2(0); %steering compliance for steering gear (back) [rad/(ft lb)]

%drag parameters; 
%do not use derivatives DLADV, DYADV, DNADV, and consequently AEROVEL
p.rho = 1.204; %air density at 1atm and 20 degree celsius [kg/m^3]
p.A = ft2_IN_m2(18.0234); %frontal area of vehicle used for longitudinal drag [m^2]
p.C_d = 0.5; %longitudinal drag coefficient [-]
p.K_tl = ft_IN_m(1); %tire drag expressed in rolling distance [m]

%suspension parameters
p.K_sf = lbs_ft_IN_N_m(1675); %suspension spring rate (front) [N/m] 
p.K_sdf = lb_sec_ft_IN_N_s_m(122.3966); %suspension damping rate (front) [N s/m] 
p.K_sr = lbs_ft_IN_N_m(1345); %suspension spring rate (rear) [N/m] 
p.K_sdr = lb_sec_ft_IN_N_s_m(112.9981); %suspension damping rate (rear) [N s/m] 

%geometric parameters
p.T_f = ft_IN_m(4.55); %track width front [m]
p.T_r = ft_IN_m(4.475); %track width rear [m]
p.h_cg = ft_IN_m(1.886053); %center of gravity height of total mass [m]

p.K_bs = lbs_ft_IN_N_m(3350); %bump stop spring rate equivalent at each wheel [N/m]
p.h_bs = ft_IN_m(0.25); %equivalent suspension clearance at bump stop at each wheel [m]
p.K_tsf = ft_lb_rad_IN_N_m_rad(-5100.155); %auxiliary torsion roll stiffness per axle (normally negative) (front) [N m/rad]
p.K_tsr = ft_lb_rad_IN_N_m_rad(-1949.82); %auxiliary torsion roll stiffness per axle (normally negative) (rear) [N m/rad]
p.K_ras = lbs_ft_IN_N_m(12000); %lateral spring rate at compliant compliant pin joint between M_s and M_u [N/m]
p.K_rad = lb_sec_ft_IN_N_s_m(700); % damping rate at compliant compliant pin joint between M_s and M_u [N s/m] 
p.K_zt = lbs_ft_IN_N_m(10842.89); % vertical spring rate of tire [N/m]

p.h_raf = ft_IN_m(0); %height of roll axis above ground (front) [m]
p.h_rar = ft_IN_m(0); %height of roll axis above ground (rear) [m]

p.h_s = ft_IN_m(2.01355); %M_s center of gravity above ground [m]

p.I_uf = lb_ft_sec2_IN_kg_m2(22.62345); %moment of inertia for unsprung mass about x-axis (front) [kg m^2]
p.I_ur = lb_ft_sec2_IN_kg_m2(21.88377); %moment of inertia for unsprung mass about x-axis (rear) [kg m^2]
p.I_y_w = []; %wheel inertia, to be determined [kg m^2]

p.K_lt = ft_lb_IN_m_N(2.397884e-4); %lateral compliance rate of tire, wheel, and suspension, per tire [m/N]

p.R_w = 0.344; %effective wheel/tire radius; chosen as tire rolling radius RR; taken from ADAMS documentation [m]
%p.R_w = ft_IN_m(0.921); %effective wheel/tire radius; chosen as tire rolling radius RR [m]
p.I_y_w = 1.7; %inertia of tire from internet forum for 235/65 R 17 [kg m^2]

p.xacc = [];
p.zacc = [];

p.dragc = [];
p.lens = [];
p.lm = [];
p.kbtf = [];
p.kvb = [];
p.kmb = [];
p.fbpvl = [];
p.xi_s = 0.5; %damping ratio for steering system lag [-]
p.omega_n_s = 70; %natural frequency for second order steering system lag [rad/sec]
p.kcf = [];
p.lso = [];
p.K_lagv = 16.5; %tire side force lag modifier for low speed operation [rad/sec]

%suspension parameters
p.D_f = rad_ft_IN_rad_m(-0.12); %[rad/m]
p.D_r = rad_ft_IN_rad_m(-0.276); %[rad/m]
p.E_f = 0; %[needs conversion if nonzero]
p.E_r = 0; %[needs conversion if nonzero]
p.K_slf = 0.105; %[-]
p.K_slr = 0.22; %[-]
p.L_saf = ft_IN_m(0.95); %[m]
p.L_sar = ft_IN_m(3.44); %[m]
p.K_sadf = 0; %[-]
p.K_sadr = 0.5; %[-]


%tire parameters from ADAMS handbook
%longitudinal coefficients
p.tire.p_cx1 = 1.6411; %Shape factor Cfx for longitudinal force
p.tire.p_dx1 = 1.1739; %Longitudinal friction Mux at Fznom
p.tire.p_dx2 = -0.16395; %Variation of friction Mux with load
p.tire.p_dx3 = 0; %Variation of friction Mux with camber
p.tire.p_ex1 = 0.46403; %Longitudinal curvature Efx at Fznom
p.tire.p_ex2 = 0.25022; %Variation of curvature Efx with load
p.tire.p_ex3 = 0.067842; %Variation of curvature Efx with load squared
p.tire.p_ex4 = -3.7604e-005; %Factor in curvature Efx while driving
p.tire.p_kx1 = 22.303; %Longitudinal slip stiffness Kfx/Fz at Fznom
p.tire.p_kx2 = 0.48896; %Variation of slip stiffness Kfx/Fz with load
p.tire.p_kx3 = 0.21253; %Exponent in slip stiffness Kfx/Fz with load
p.tire.p_hx1 = 0.0012297; %Horizontal shift Shx at Fznom
p.tire.p_hx2 = 0.0004318; %Variation of shift Shx with load
p.tire.p_vx1 = -8.8098e-006; %Vertical shift Svx/Fz at Fznom
p.tire.p_vx2 = 1.862e-005; %Variation of shift Svx/Fz with load
p.tire.r_bx1 = 13.276; %Slope factor for combined slip Fx reduction
p.tire.r_bx2 = -13.778; %Variation of slope Fx reduction with kappa
p.tire.r_cx1 = 1.2568; %Shape factor for combined slip Fx reduction
p.tire.r_ex1 = 0.65225; %Curvature factor of combined Fx
p.tire.r_ex2 = -0.24948; %Curvature factor of combined Fx with load
p.tire.r_hx1 = 0.0050722; %Shift factor for combined slip Fx reduction
p.tire.p_tx1 = 2.3657; %Relaxation length SigKap0/Fz at Fznom
p.tire.p_tx2 = 1.4112; %Variation of SigKap0/Fz with load
p.tire.p_tx3 = 0.56626; %Variation of SigKap0/Fz with exponent of load

%lateral coefficients
p.tire.p_cy1 = 1.3507; %Shape factor Cfy for lateral forces
p.tire.p_dy1 = 1.0489; %Lateral friction Muy
p.tire.p_dy2 = -0.18033; %Variation of friction Muy with load
p.tire.p_dy3 = -2.8821; %Variation of friction Muy with squared camber
p.tire.p_ey1 = -0.0074722; %Lateral curvature Efy at Fznom
p.tire.p_ey2 = -0.0063208; %Variation of curvature Efy with load
p.tire.p_ey3 = -9.9935; %Zero order camber dependency of curvature Efy
p.tire.p_ey4 = -760.14; %Variation of curvature Efy with camber
p.tire.p_ky1 = -21.92; %Maximum value of stiffness Kfy/Fznom
p.tire.p_ky2 = 2.0012; %Load at which Kfy reaches maximum value
p.tire.p_ky3 = -0.024778; %Variation of Kfy/Fznom with camber
p.tire.p_hy1 = 0.0026747; %Horizontal shift Shy at Fznom
p.tire.p_hy2 = 8.9094e-005; %Variation of shift Shy with load
p.tire.p_hy3 = 0.031415; %Variation of shift Shy with camber
p.tire.p_vy1 = 0.037318; %Vertical shift in Svy/Fz at Fznom
p.tire.p_vy2 = -0.010049; %Variation of shift Svy/Fz with load
p.tire.p_vy3 = -0.32931; %Variation of shift Svy/Fz with camber
p.tire.p_vy4 = -0.69553; %Variation of shift Svy/Fz with camber and load
p.tire.r_by1 = 7.1433; %Slope factor for combined Fy reduction
p.tire.r_by2 = 9.1916; %Variation of slope Fy reduction with alpha
p.tire.r_by3 = -0.027856; %Shift term for alpha in slope Fy reduction
p.tire.r_cy1 = 1.0719; %Shape factor for combined Fy reduction
p.tire.r_ey1 = -0.27572; %Curvature factor of combined Fy
p.tire.r_ey2 = 0.32802; %Curvature factor of combined Fy with load
p.tire.r_hy1 = 5.7448e-006; %Shift factor for combined Fy reduction
p.tire.r_hy2 = -3.1368e-005; %Shift factor for combined Fy reduction with load
p.tire.r_vy1 = -0.027825; %Kappa induced side force Svyk/Muy*Fz at Fznom
p.tire.r_vy2 = 0.053604; %Variation of Svyk/Muy*Fz with load
p.tire.r_vy3 = -0.27568; %Variation of Svyk/Muy*Fz with camber
p.tire.r_vy4 = 12.12; %Variation of Svyk/Muy*Fz with alpha
p.tire.r_vy5 = 1.9; %Variation of Svyk/Muy*Fz with kappa
p.tire.r_vy6 = -10.704; %Variation of Svyk/Muy*Fz with atan(kappa)
p.tire.p_ty1 = 2.1439; %Peak value of relaxation length SigAlp0/R0
p.tire.p_ty2 = 1.9829; %Value of Fz/Fznom where SigAlp0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%load parameters
g = 9.81; %[m/s^2]

%create equivalent bicycle parameters
%mu = p.tire.p_dy1;
mu = 0.9; %<-- hard coded
C_Sf = -p.tire.p_ky1/p.tire.p_dy1; %<--corrected
C_Sr = -p.tire.p_ky1/p.tire.p_dy1; %<--corrected
lf = p.a;
lr = p.b;
h = p.h_s;
%h = 0;
m = p.m;
I = p.I_z;

%states
%x1 = β slip angle at vehicle center
%x2 = Ψ yaw angle
%x3 = Ψ yaw rate
%x4 = u velocity in x-direction
%x5 = s_x x-position in a global coordinate system
%x6 = s_y y-position in a global coordinate system

%u(1) = delta_w steering angle of front wheels
%u(2) = ax longitudinal acceleration



%system dynamics
f(1,1) = (mu/(x(4)^2*(lr+lf))*(C_Sr*(g*lf + (x(8)+u(2))*h)*lr - C_Sf*(g*lr - (x(8)+u(2))*h)*lf)-1)*x(3) ...
    -mu/(x(4)*(lr+lf))*(C_Sr*(g*lf + (x(8)+u(2))*h) + C_Sf*(g*lr-(x(8)+u(2))*h))*x(1) ...
    +mu/(x(4)*(lr+lf))*(C_Sf*(g*lr-(x(8)+u(2))*h))*(x(7)+u(1));
f(2,1) = x(3);
f(3,1) = -mu*m/(x(4)*I*(lr+lf))*(lf^2*C_Sf*(g*lr-(x(8)+u(2))*h) + lr^2*C_Sr*(g*lf + (x(8)+u(2))*h))*x(3) ...
    +mu*m/(I*(lr+lf))*(lr*C_Sr*(g*lf + (x(8)+u(2))*h) - lf*C_Sf*(g*lr - (x(8)+u(2))*h))*x(1) ...
    +mu*m/(I*(lr+lf))*lf*C_Sf*(g*lr - (x(8)+u(2))*h)*(x(7)+u(1));
f(4,1) = (x(8)+u(2));
f(5,1) = x(4)*cos(x(1) + x(2));
f(6,1) = x(4)*sin(x(1) + x(2));
f(7,1) = 0;
f(8,1) = 0;
f(9,1) = 0;
f(10,1) = 0;
f(11,1) = 0;
f(12,1) = 0;
f(13,1) = 0;
f(14,1) = 0;
f(15,1) = 0;
f(16,1) = 0;
f(17,1) = 0;
f(18,1) = 0;





%system dynamics
% f(1,1) = (mu/(x(4)^2*(lr+lf))*(C_Sr*(g*lf + u(2)*h)*lr - C_Sf*(g*lr - u(2)*h)*lf)-1)*x(3) ...
%     -mu/(x(4)*(lr+lf))*(C_Sr*(g*lf + u(2)*h) + C_Sf*(g*lr-u(2)*h))*x(1) ...
%     +mu/(x(4)*(lr+lf))*(C_Sf*(g*lr-u(2)*h))*u(1);
% f(2,1) = x(3);
% f(3,1) = -mu*m/(x(4)*I*(lr+lf))*(lf^2*C_Sf*(g*lr-u(2)*h) + lr^2*C_Sr*(g*lf + u(2)*h))*x(3) ...
%     +mu*m/(I*(lr+lf))*(lr*C_Sr*(g*lf + u(2)*h) - lf*C_Sf*(g*lr - u(2)*h))*x(1) ...
%     +mu*m/(I*(lr+lf))*lf*C_Sf*(g*lr - u(2)*h)*u(1);
% f(4,1) = u(2);
% f(5,1) = x(4)*cos(x(1) + x(2));
% f(6,1) = x(4)*sin(x(1) + x(2));

%------------- END OF CODE --------------





function post_val = lb_sec2_ft_IN_kg(prev_val)

% 1lb is 4.4482216152605 N
% 1ft is 0.3048 m

post_val = 4.4482216152605/0.3048*prev_val;


function post_val = ft_IN_m(prev_val)
%original: [ft]
%new: [m]

% 1ft is 0.3048 m

post_val = 0.3048*prev_val;


function post_val = lb_ft_sec2_IN_kg_m2(prev_val)

%[kg m^2] = [N m sec^2]

% 1lb is 4.4482216152605 N
% 1ft is 0.3048 m

post_val = 4.4482216152605*0.3048*prev_val;

function post_val = rad_ft_lb_IN_rad_sec2_kg_m2(prev_val)

%original: [rad/(ft lb)]
%new: [rad/(N m)] = [rad s^2/(kg m^2)]

% 1lb is 4.4482216152605 N
% 1ft is 0.3048 m

post_val = 1/(4.4482216152605*0.3048)*prev_val;

function post_val = ft2_IN_m2(prev_val)
%original: [ft^2]
%new: [m^2]

% 1ft is 0.3048 m

post_val = 0.3048^2*prev_val;

function post_val = lbs_ft_IN_N_m(prev_val)
%original: [lbs/ft]
%new: [N/m]

% 1lbs is 0.45359237 kg
% 1kg is around 9.81 N assuming being close to sea level
% 1ft is 0.3048 m

post_val = 0.45359237*9.81/0.3048*prev_val;

function post_val = lb_sec_ft_IN_N_s_m(prev_val)
%original: [lb sec/ft]
%new: [N sec/m]

% 1lb is 4.4482216152605 N
% 1ft is 0.3048 m

post_val = 4.4482216152605/0.3048*prev_val;

function post_val = ft_lb_rad_IN_N_m_rad(prev_val)
%original: [lb ft/rad]
%new: [N m/rad]

% 1lb is 4.4482216152605 N
% 1ft is 0.3048 m

post_val = 4.4482216152605*0.3048*prev_val;

function post_val = ft_lb_IN_m_N(prev_val)
%original: [ft/lb]
%new: [m/N]

% 1lb is 4.4482216152605 N
% 1ft is 0.3048 m

post_val = 0.3048/4.4482216152605*prev_val;

function post_val = rad_ft_IN_rad_m(prev_val)
%original: [rad/ft]
%new: [rad/m]

% 1ft is 0.3048 m

post_val = 1/0.3048*prev_val;
