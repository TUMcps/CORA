function p = BMWparameters()
% BMWparameters - parameter for a BMW vehicle (see [1])
%
% Syntax:  
%    p = BMWparameters()
%
% Inputs:
%    ---
%
% Outputs:
%    p - struct storing the parameters
%
% References:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: DOTBicycleDynamics_controlled_BMW

% Author:       Matthias Althoff
% Written:      16-August-2011
% Last update:  16-June-2023
% Last revision:---

%------------- BEGIN CODE --------------

%friction
p.mu = 1;

%vehicle body dimensions
p.l = 4.508; %vehicle length [m] (with US bumpers)
p.w = 1.610; %vehicle width [m]

%steering constraints
p.steering.min = -1.066; %minimum steering angle [rad]
p.steering.max = 1.066; %maximum steering angle [rad]
p.steering.v_min = -0.4; %minimum steering velocity [rad/s]
p.steering.v_max = 0.4; %maximum steering velocity [rad/s]

%longitudinal constraints
p.longitudinal.v_min = -13.6; %minimum velocity [m/s]
p.longitudinal.v_max = 50.8; %minimum velocity [m/s]
p.longitudinal.v_switch = 7.319; %switching velocity [m/s]
p.longitudinal.a_max = 11.5; %maximum absolute acceleration [m/s^2]

%masses
p.m = lb_sec2_ft_IN_kg(74.91452); %vehicle mass [kg]; MASS
p.m_s = lb_sec2_ft_IN_kg(66.17221); %sprung mass [kg]; SMASS
p.m_uf = lb_sec2_ft_IN_kg(4.371153); %unsprung mass front [kg]; UMASSF
p.m_ur = lb_sec2_ft_IN_kg(4.371153); %unsprung mass rear [kg]; UMASSR

%axes distances
p.a = ft_IN_m(3.793293); %distance from spring mass center of gravity to front axle [m]; LENA
p.b = ft_IN_m(4.667707); %distance from spring mass center of gravity to rear axle [m]; LENB

%moments of inertia of sprung mass
p.I_Phi_s = lb_ft_sec2_IN_kg_m2(152.871); %moment of inertia for sprung mass in roll [kg m^2]; IXS
p.I_y_s = lb_ft_sec2_IN_kg_m2(1154.888); %moment of inertia for sprung mass in pitch [kg m^2]; IYS
p.I_z = lb_ft_sec2_IN_kg_m2(1321.416); %moment of inertia for sprung mass in yaw [kg m^2]; IZZ
p.I_xz_s = lb_ft_sec2_IN_kg_m2(0); %moment of inertia cross product [kg m^2]; IXZ

%suspension parameters
p.K_sf = lbs_ft_IN_N_m(1675); %suspension spring rate (front) [N/m]; KSF
p.K_sdf = lb_sec_ft_IN_N_s_m(122.3966); %suspension damping rate (front) [N s/m]; KSDF
p.K_sr = lbs_ft_IN_N_m(1345); %suspension spring rate (rear) [N/m]; KSR
p.K_sdr = lb_sec_ft_IN_N_s_m(112.9981); %suspension damping rate (rear) [N s/m]; KSDR 

%geometric parameters
p.T_f = ft_IN_m(4.55); %track width front [m]; TRWF
p.T_r = ft_IN_m(4.475); %track width rear [m]; TRWB
p.K_ras = lbs_ft_IN_N_m(12000); %lateral spring rate at compliant compliant pin joint between M_s and M_u [N/m]; KRAS

p.K_tsf = ft_lb_rad_IN_N_m_rad(-5100.155); %auxiliary torsion roll stiffness per axle (normally negative) (front) [N m/rad]; KTSF
p.K_tsr = ft_lb_rad_IN_N_m_rad(-1949.82); %auxiliary torsion roll stiffness per axle (normally negative) (rear) [N m/rad]; KTSR
p.K_rad = lb_sec_ft_IN_N_s_m(700); % damping rate at compliant compliant pin joint between M_s and M_u [N s/m]; KRADP
p.K_zt = lbs_ft_IN_N_m(10842.89); % vertical spring rate of tire [N/m]; TSPRINGR

p.h_cg = ft_IN_m(1.886053); %center of gravity height of total mass [m]; HCG (mainly required for conversion to other vehicle models)
p.h_raf = ft_IN_m(0); %height of roll axis above ground (front) [m]; HRAF
p.h_rar = ft_IN_m(0); %height of roll axis above ground (rear) [m]; HRAR

p.h_s = ft_IN_m(2.01355); %M_s center of gravity above ground [m]; HS

p.I_uf = lb_ft_sec2_IN_kg_m2(22.62345); %moment of inertia for unsprung mass about x-axis (front) [kg m^2]; IXUF
p.I_ur = lb_ft_sec2_IN_kg_m2(21.88377); %moment of inertia for unsprung mass about x-axis (rear) [kg m^2]; IXUR
p.I_y_w = 1.7; %wheel inertia, from internet forum for 235/65 R 17 [kg m^2]

p.K_lt = ft_lb_IN_m_N(2.397884e-4); %lateral compliance rate of tire, wheel, and suspension, per tire [m/N]; KLT
p.R_w = 0.344; %effective wheel/tire radius; chosen as tire rolling radius RR; taken from ADAMS documentation [m]

%split of brake and engine torque
p.T_sb = 0.66;
p.T_se = 0;

%suspension parameters
p.D_f = rad_ft_IN_rad_m(-0.12); %[rad/m]; DF
p.D_r = rad_ft_IN_rad_m(-0.276); %[rad/m]; DR
p.E_f = 0; %[needs conversion if nonzero]; EF
p.E_r = 0; %[needs conversion if nonzero]; ER


%tire parameters from ADAMS handbook
%longitudinal coefficients
p.tire.p_cx1 = 1.6411; %Shape factor Cfx for longitudinal force
p.tire.p_dx1 = 1.1739; %Longitudinal friction Mux at Fznom
p.tire.p_dx3 = 0; %Variation of friction Mux with camber
p.tire.p_ex1 = 0.46403; %Longitudinal curvature Efx at Fznom
p.tire.p_kx1 = 22.303; %Longitudinal slip stiffness Kfx/Fz at Fznom
p.tire.p_hx1 = 0.0012297; %Horizontal shift Shx at Fznom
p.tire.p_vx1 = -8.8098e-006; %Vertical shift Svx/Fz at Fznom
p.tire.r_bx1 = 13.276; %Slope factor for combined slip Fx reduction
p.tire.r_bx2 = -13.778; %Variation of slope Fx reduction with kappa
p.tire.r_cx1 = 1.2568; %Shape factor for combined slip Fx reduction
p.tire.r_ex1 = 0.65225; %Curvature factor of combined Fx
p.tire.r_hx1 = 0.0050722; %Shift factor for combined slip Fx reduction

%lateral coefficients
p.tire.p_cy1 = 1.3507; %Shape factor Cfy for lateral forces
p.tire.p_dy1 = 1.0489; %Lateral friction Muy
p.tire.p_dy3 = -2.8821; %Variation of friction Muy with squared camber
p.tire.p_ey1 = -0.0074722; %Lateral curvature Efy at Fznom
p.tire.p_ky1 = -21.92; %Maximum value of stiffness Kfy/Fznom
p.tire.p_hy1 = 0.0026747; %Horizontal shift Shy at Fznom
p.tire.p_hy3 = 0.031415; %Variation of shift Shy with camber
p.tire.p_vy1 = 0.037318; %Vertical shift in Svy/Fz at Fznom
p.tire.p_vy3 = -0.32931; %Variation of shift Svy/Fz with camber
p.tire.r_by1 = 7.1433; %Slope factor for combined Fy reduction
p.tire.r_by2 = 9.1916; %Variation of slope Fy reduction with alpha
p.tire.r_by3 = -0.027856; %Shift term for alpha in slope Fy reduction
p.tire.r_cy1 = 1.0719; %Shape factor for combined Fy reduction
p.tire.r_ey1 = -0.27572; %Curvature factor of combined Fy
p.tire.r_hy1 = 5.7448e-006; %Shift factor for combined Fy reduction
p.tire.r_vy1 = -0.027825; %Kappa induced side force Svyk/Muy*Fz at Fznom
p.tire.r_vy3 = -0.27568; %Variation of Svyk/Muy*Fz with camber
p.tire.r_vy4 = 12.12; %Variation of Svyk/Muy*Fz with alpha
p.tire.r_vy5 = 1.9; %Variation of Svyk/Muy*Fz with kappa
p.tire.r_vy6 = -10.704; %Variation of Svyk/Muy*Fz with atan(kappa)



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