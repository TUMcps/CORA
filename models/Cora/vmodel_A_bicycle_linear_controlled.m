function dx = vmodel_A_bicycle_linear_controlled(t,x,u)
%Bicycle Model with 
% - normal force equilibrium for pitching-moments
% - linear tyre model
% state x=[X,Y,psi,vx,vy,omega]
% input u=[delta,omega_f,omega_r]

%parameters: get parameters from p vector
    %body
    m = 1750;
    J = 2500;
    L = 2.7;
    l_F = 1.43;
    l_R = L-l_F;
    h = 0.5;

    %street
    mu0 = 1;
    g = 9.81;
    
    %tires
    %B=p(7),C=p(8) - paceijca parameter in ff_abs = sin(C*atan(B*sf_abs/mu0))
    C_F = 10.4*1.3;
    C_R = 21.4*1.1;

%state
    %position
    X = x(1); %#ok<NASGU>
    Y = x(2); %#ok<NASGU>
    psi = x(3);
    
    %velocity
    vx = x(4);
    vy = x(5);
    omega = x(6);
    
    %acceleration
    Fb = x(7)*m;
    delta = x(8);
    
    
%control action
% inputs are values of the state feedback matrix R, the reference state Xn,
% and the feedforward value W

R = [u(1) u(2) u(3) u(4) u(5) u(6) u(7) u(8);...
     u(9) u(10) u(11) u(12) u(13) u(14) u(15) u(16)];
Xn = [u(17); u(18); u(19); u(20); u(21); u(22); u(23); u(24)];
W = [u(25); u(26)];
v = -R*(x-Xn)+W;
        
%calculate normal forces
Fzf = (l_R*m*g-h*Fb)/(l_R+l_F);
Fzr = m*g - Fzf;

%side-slip
sf =(vy+l_F*omega)/vx-delta;
sr =(vy-l_R*omega)/vx;

%forces
Fyf = -C_F * Fzf * sf;
Fyr = -C_R * Fzr * sr;


%ACCELERATIONS
dvx = Fb/m + vy * omega;
dvy = (Fyf+Fyr)/m - vx * omega;
domega = (l_F * Fyf - l_R * Fyr)/J;


%position
    cp = cos(psi);sp = sin(psi);
    dx(1,1) = cp * vx - sp * vy;
    dx(2,1) = sp * vx + cp * vy;
    dx(3,1) = omega;
%velocity
    dx(4,1) = dvx;
    dx(5,1) = dvy;
    dx(6,1) = domega;
%accelerationm
    dx(7,1) = v(1);
    dx(8,1) = v(2);
