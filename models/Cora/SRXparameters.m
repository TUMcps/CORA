function p = SRXparameters()
% updated: 14-March-2012, MA


%masses
p.m = 2273.4049; %vehicle mass [kg]

%axes distances
p.a = 1.2916; %distance from spring mass center of gravity to front axle [m]
p.b = 1.5151; %distance from spring mass center of gravity to rear axle [m]

%moments of inertia of sprung mass
p.I_z = 4423.0301; %moment of inertia for sprung mass in yaw [kg m^2]


%tire parameters from ADAMS handbook
%lateral coefficients
p.tire.p_dy1 = 1; %Lateral friction Muy
p.tire.p_ky1 = -107555.5915/(p.m*9.81/4); % old value: 21.92; %Maximum value of stiffness Kfy/Fznom





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