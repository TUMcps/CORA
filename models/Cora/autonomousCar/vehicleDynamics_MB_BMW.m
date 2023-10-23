function f = vehicleDynamics_MB_BMW(x,u)
% vehicleDynamics_MB_BMW - multi-body vehicle dynamics of a BMW
% reference point: center of mass
%
% Syntax:  
%    f = vehicleDynamics_MB_BMW(t,x,u,p)
%
% Inputs:
%    x - vehicle state vector
%    u - vehicle input vector
%    p - vehicle parameter structure
%
% Outputs:
%    f - right-hand side of differential equations
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      05-January-2017
% Last update:  15-December-2017
%               21-June-2023
% Last revision:---

%------------- BEGIN CODE --------------

% set gravity constant
g = 9.81; %[m/s^2]

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

% steering velocity constraints and acceleation constraints are not
% considered for now

%compute slip angle at cg
%switch to kinematic model for small velocities
if abs(x(4)) < 0.1
    beta = 0;
else
    beta = atan(x(11)/x(4));
end
vel = sqrt(x(4)^2 + x(11)^2);

%vertical tire forces
F_z_LF = (x(17) + p.R_w*(cos(x(14)) - 1) - 0.5*p.T_f*sin(x(14)))*p.K_zt;
F_z_RF = (x(17) + p.R_w*(cos(x(14)) - 1) + 0.5*p.T_f*sin(x(14)))*p.K_zt;
F_z_LR = (x(22) + p.R_w*(cos(x(19)) - 1) - 0.5*p.T_r*sin(x(19)))*p.K_zt;
F_z_RR = (x(22) + p.R_w*(cos(x(19)) - 1) + 0.5*p.T_r*sin(x(19)))*p.K_zt;

%obtain individual tire speeds
u_w_lf = (x(4) + 0.5*p.T_f*x(6))*cos(x(3)) + (x(11) + p.a*x(6))*sin(x(3));
u_w_rf = (x(4) - 0.5*p.T_f*x(6))*cos(x(3)) + (x(11) + p.a*x(6))*sin(x(3));
u_w_lr = x(4) + 0.5*p.T_r*x(6);
u_w_rr = x(4) - 0.5*p.T_r*x(6);

%compute longitudinal slip
%switch to kinematic model for small velocities
if abs(x(4)) < 0.1
    s_lf = 0;
    s_rf = 0;
    s_lr = 0;
    s_rr = 0;    
else
    s_lf = 1 - p.R_w*x(24)/u_w_lf;
    s_rf = 1 - p.R_w*x(25)/u_w_rf;
    s_lr = 1 - p.R_w*x(26)/u_w_lr;
    s_rr = 1 - p.R_w*x(27)/u_w_rr;
end

%lateral slip angles
%switch to kinematic model for small velocities
if abs(x(4)) < 0.1
    alpha_LF = 0;
    alpha_RF = 0;
    alpha_LR = 0;
    alpha_RR = 0;
else
    alpha_LF = atan((x(11) + p.a*x(6) - x(15)*(p.R_w - x(17)))/(x(4) + 0.5*p.T_f*x(6))) - x(3);
    alpha_RF = atan((x(11) + p.a*x(6) - x(15)*(p.R_w - x(17)))/(x(4) - 0.5*p.T_f*x(6))) - x(3);
    alpha_LR = atan((x(11) - p.b*x(6) - x(20)*(p.R_w - x(22)))/(x(4) + 0.5*p.T_r*x(6)));
    alpha_RR = atan((x(11) - p.b*x(6) - x(20)*(p.R_w - x(22)))/(x(4) - 0.5*p.T_r*x(6)));
end

%auxiliary suspension movement
z_SLF = (p.h_s - p.R_w + x(17) - x(12))/cos(x(7)) - p.h_s + p.R_w + p.a*x(9) + 0.5*(x(7) - x(14))*p.T_f;
z_SRF = (p.h_s - p.R_w + x(17) - x(12))/cos(x(7)) - p.h_s + p.R_w + p.a*x(9) - 0.5*(x(7) - x(14))*p.T_f;
z_SLR = (p.h_s - p.R_w + x(22) - x(12))/cos(x(7)) - p.h_s + p.R_w - p.b*x(9) + 0.5*(x(7) - x(19))*p.T_r;
z_SRR = (p.h_s - p.R_w + x(22) - x(12))/cos(x(7)) - p.h_s + p.R_w - p.b*x(9) - 0.5*(x(7) - x(19))*p.T_r;

dz_SLF = x(18) - x(13) + p.a*x(10) + 0.5*(x(8) - x(15))*p.T_f;
dz_SRF = x(18) - x(13) + p.a*x(10) - 0.5*(x(8) - x(15))*p.T_f;
dz_SLR = x(23) - x(13) - p.b*x(10) + 0.5*(x(8) - x(20))*p.T_r;
dz_SRR = x(23) - x(13) - p.b*x(10) - 0.5*(x(8) - x(20))*p.T_r;

%camber angles
gamma_LF = x(7) + p.D_f*z_SLF + p.E_f*(z_SLF)^2;
gamma_RF = x(7) - p.D_f*z_SRF - p.E_f*(z_SRF)^2;
gamma_LR = x(7) + p.D_r*z_SLR + p.E_r*(z_SLR)^2;
gamma_RR = x(7) - p.D_r*z_SRR - p.E_r*(z_SRR)^2;

%compute longitudinal tire forces using the magic formula for pure slip
F0_x_LF = mFormulaLongitudinal(s_lf, gamma_LF, F_z_LF, p.tire);
F0_x_RF = mFormulaLongitudinal(s_rf, gamma_RF, F_z_RF, p.tire);
F0_x_LR = mFormulaLongitudinal(s_lr, gamma_LR, F_z_LR, p.tire);
F0_x_RR = mFormulaLongitudinal(s_rr, gamma_RR, F_z_RR, p.tire);

%compute lateral tire forces using the magic formula for pure slip
[F0_y_LF, mu_y_LF] = mFormulaLateral(alpha_LF, gamma_LF, F_z_LF, p.tire);
[F0_y_RF, mu_y_RF] = mFormulaLateral(alpha_RF, gamma_RF, F_z_RF, p.tire);
[F0_y_LR, mu_y_LR] = mFormulaLateral(alpha_LR, gamma_LR, F_z_LR, p.tire);
[F0_y_RR, mu_y_RR] = mFormulaLateral(alpha_RR, gamma_RR, F_z_RR, p.tire);

%compute longitudinal tire forces using the magic formula for combined slip
F_x_LF = mFormulaLongitudinalComb(s_lf, alpha_LF, F0_x_LF, p.tire);
F_x_RF = mFormulaLongitudinalComb(s_rf, alpha_RF, F0_x_RF, p.tire);
F_x_LR = mFormulaLongitudinalComb(s_lr, alpha_LR, F0_x_LR, p.tire);
F_x_RR = mFormulaLongitudinalComb(s_rr, alpha_RR, F0_x_RR, p.tire);

%compute lateral tire forces using the magic formula for combined slip
F_y_LF = mFormulaLateralComb(s_lf, alpha_LF, gamma_LF, mu_y_LF, F_z_LF, F0_y_LF, p.tire);
F_y_RF = mFormulaLateralComb(s_rf, alpha_RF, gamma_RF, mu_y_RF, F_z_RF, F0_y_RF, p.tire);
F_y_LR = mFormulaLateralComb(s_lr, alpha_LR, gamma_LR, mu_y_LR, F_z_LR, F0_y_LR, p.tire);
F_y_RR = mFormulaLateralComb(s_rr, alpha_RR, gamma_RR, mu_y_RR, F_z_RR, F0_y_RR, p.tire);

%auxiliary movements for compliant joint equations
delta_z_f = p.h_s - p.R_w + x(17) - x(12);
delta_z_r = p.h_s - p.R_w + x(22) - x(12);

delta_phi_f = x(7) - x(14);
delta_phi_r = x(7) - x(19);

dot_delta_phi_f = x(8) - x(15);
dot_delta_phi_r = x(8) - x(20);

dot_delta_z_f = x(18) - x(13);
dot_delta_z_r = x(23) - x(13);

dot_delta_y_f = x(11) + p.a*x(6) - x(16);
dot_delta_y_r = x(11) - p.b*x(6) - x(21);

delta_f = delta_z_f*sin(x(7)) - x(28)*cos(x(7)) - (p.h_raf - p.R_w)*sin(delta_phi_f);
delta_r = delta_z_r*sin(x(7)) - x(29)*cos(x(7)) - (p.h_rar - p.R_w)*sin(delta_phi_r);

dot_delta_f = (delta_z_f*cos(x(7)) + x(28)*sin(x(7)))*x(8)...
              + dot_delta_z_f*sin(x(7)) - dot_delta_y_f*cos(x(7)) - (p.h_raf - p.R_w)*cos(delta_phi_f)*dot_delta_phi_f;
dot_delta_r = (delta_z_r*cos(x(7)) + x(29)*sin(x(7)))*x(8)...
              + dot_delta_z_r*sin(x(7)) - dot_delta_y_r*cos(x(7)) - (p.h_rar - p.R_w)*cos(delta_phi_r)*dot_delta_phi_r;    
          
%compliant joint forces
F_RAF = delta_f*p.K_ras + dot_delta_f*p.K_rad;
F_RAR = delta_r*p.K_ras + dot_delta_r*p.K_rad;

%auxiliary suspension forces (bump stop neglected; squat/lift forces neglected)
F_SLF = p.m_s*g*p.b/(2*(p.a+p.b)) - z_SLF*p.K_sf - dz_SLF*p.K_sdf + (x(7) - x(14))*p.K_tsf/p.T_f;

F_SRF = p.m_s*g*p.b/(2*(p.a+p.b)) - z_SRF*p.K_sf - dz_SRF*p.K_sdf - (x(7) - x(14))*p.K_tsf/p.T_f;

F_SLR = p.m_s*g*p.a/(2*(p.a+p.b)) - z_SLR*p.K_sr - dz_SLR*p.K_sdr + (x(7) - x(19))*p.K_tsr/p.T_r;

F_SRR = p.m_s*g*p.a/(2*(p.a+p.b)) - z_SRR*p.K_sr - dz_SRR*p.K_sdr - (x(7) - x(19))*p.K_tsr/p.T_r;


%auxiliary variables sprung mass
sumX = F_x_LR + F_x_RR + (F_x_LF + F_x_RF)*cos(x(3)) - (F_y_LF + F_y_RF)*sin(x(3));

sumN = (F_y_LF + F_y_RF)*p.a*cos(x(3)) + (F_x_LF + F_x_RF)*p.a*sin(x(3))...
       + (F_y_RF - F_y_LF)*0.5*p.T_f*sin(x(3)) + (F_x_LF - F_x_RF)*0.5*p.T_f*cos(x(3))...
       + (F_x_LR - F_x_RR)*0.5*p.T_r - (F_y_LR + F_y_RR)*p.b;
   
sumY_s = (F_RAF + F_RAR)*cos(x(7)) + (F_SLF + F_SLR + F_SRF + F_SRR)*sin(x(7));

sumL = 0.5*F_SLF*p.T_f + 0.5*F_SLR*p.T_r - 0.5*F_SRF*p.T_f - 0.5*F_SRR*p.T_r...
       - F_RAF/cos(x(7))*(p.h_s - x(12) - p.R_w + x(17) - (p.h_raf - p.R_w)*cos(x(14)))...
       - F_RAR/cos(x(7))*(p.h_s - x(12) - p.R_w + x(22) - (p.h_rar - p.R_w)*cos(x(19)));
   
sumZ_s = (F_SLF + F_SLR + F_SRF + F_SRR)*cos(x(7)) - (F_RAF + F_RAR)*sin(x(7));

sumM_s = p.a*(F_SLF + F_SRF) - p.b*(F_SLR + F_SRR) + ((F_x_LF + F_x_RF)*cos(x(3))...
       - (F_y_LF + F_y_RF)*sin(x(3)) + F_x_LR + F_x_RR)*(p.h_s - x(12));

%auxiliary variables unsprung mass
sumL_uf = 0.5*F_SRF*p.T_f - 0.5*F_SLF*p.T_f - F_RAF*(p.h_raf - p.R_w)...
          + F_z_LF*(p.R_w*sin(x(14)) + 0.5*p.T_f*cos(x(14)) - p.K_lt*F_y_LF)...
          - F_z_RF*(-p.R_w*sin(x(14)) + 0.5*p.T_f*cos(x(14)) + p.K_lt*F_y_RF)...
          - ((F_y_LF + F_y_RF)*cos(x(3)) + (F_x_LF + F_x_RF)*sin(x(3)))*(p.R_w - x(17));
      
sumL_ur = 0.5*F_SRR*p.T_r - 0.5*F_SLR*p.T_r - F_RAR*(p.h_rar - p.R_w)...
          + F_z_LR*(p.R_w*sin(x(19)) + 0.5*p.T_r*cos(x(19)) - p.K_lt*F_y_LR)...
          - F_z_RR*(-p.R_w*sin(x(19)) + 0.5*p.T_r*cos(x(19)) + p.K_lt*F_y_RR)...
          - (F_y_LR + F_y_RR)*(p.R_w - x(22));   
      
sumZ_uf = F_z_LF + F_z_RF + F_RAF*sin(x(7)) - (F_SLF + F_SRF)*cos(x(7));

sumZ_ur = F_z_LR + F_z_RR + F_RAR*sin(x(7)) - (F_SLR + F_SRR)*cos(x(7));

sumY_uf = (F_y_LF + F_y_RF)*cos(x(3)) + (F_x_LF + F_x_RF)*sin(x(3))...
          - F_RAF*cos(x(7)) - (F_SLF + F_SRF)*sin(x(7));
      
sumY_ur = (F_y_LR + F_y_RR)...
          - F_RAR*cos(x(7)) - (F_SLR + F_SRR)*sin(x(7)); 
      
      
%dynamics common with single-track model
%switch to kinematic model for small velocities
if abs(x(4)) < 0.1
    % use kinematic model with ref point at center of mass
    lwb = p.a + p.b;
    %system dynamics
    f(1:5,1) = vehicleDynamics_KS_cog(x(1:5),u,p);
    d_beta = (p.b * u(1)) / (lwb*cos(x(3))^2 * (1 + (tan(x(3))^2 * p.b/lwb)^2));
    dd_psi = 1/lwb * (u(2)*cos(x(7))*tan(x(3)) - x(3)*sin(x(7))*d_beta*tan(x(3)) + x(4)*cos(x(7))*u(1)/cos(x(3))^2);
    f(6,1) = dd_psi;
else
    f(1,1) = cos(beta + x(5))*vel;
    f(2,1) = sin(beta + x(5))*vel;
    f(3,1) = u(1);
    f(4,1) = 1/p.m*sumX + x(6)*x(11);
    f(5,1) = x(6);
    f(6,1) = 1/(p.I_z - (p.I_xz_s)^2/p.I_Phi_s)*(sumN + p.I_xz_s/p.I_Phi_s*sumL);
end


% remaining sprung mass dynamics
f(7,1) = x(8);
f(8,1) = 1/(p.I_Phi_s - (p.I_xz_s)^2/p.I_z)*(p.I_xz_s/p.I_z*sumN + sumL);
f(9,1) = x(10);
f(10,1) = 1/p.I_y_s*sumM_s;
f(11,1) = 1/p.m_s*sumY_s - x(6)*x(4);
f(12,1) = x(13);
f(13,1) = g - 1/p.m_s*sumZ_s;

%unsprung mass dynamics (front)
f(14,1) = x(15);
f(15,1) = 1/p.I_uf*sumL_uf;
f(16,1) = 1/p.m_uf*sumY_uf - x(6)*x(4);
f(17,1) = x(18);
f(18,1) = g - 1/p.m_uf*sumZ_uf;

%unsprung mass dynamics (rear)
f(19,1) = x(20);
f(20,1) = 1/p.I_ur*sumL_ur;
f(21,1) = 1/p.m_ur*sumY_ur - x(6)*x(4);
f(22,1) = x(23);
f(23,1) = g - 1/p.m_ur*sumZ_ur;

%convert acceleration input to brake and engine torque
if u(2)>0
    T_B = 0;
    T_E = p.m*p.R_w*u(2);
else
    T_B = p.m*p.R_w*u(2);
    T_E = 0;
end

%wheel dynamics (p.T; new parameter for torque splitting)
f(24,1) = 1/p.I_y_w*(-p.R_w*F_x_LF + 0.5*p.T_sb*T_B + 0.5*p.T_se*T_E);
f(25,1) = 1/p.I_y_w*(-p.R_w*F_x_RF + 0.5*p.T_sb*T_B + 0.5*p.T_se*T_E);
f(26,1) = 1/p.I_y_w*(-p.R_w*F_x_LR + 0.5*(1-p.T_sb)*T_B + 0.5*(1-p.T_se)*T_E);
f(27,1) = 1/p.I_y_w*(-p.R_w*F_x_RR + 0.5*(1-p.T_sb)*T_B + 0.5*(1-p.T_se)*T_E);

%negative wheel spin forbidden
for iState = 24:27
    if x(iState)<0
       x(iState) = 0;
       f(iState,1) = 0; 
    end
end

%compliant joint equations
f(28,1) = dot_delta_y_f;
f(29,1) = dot_delta_y_r;

end

function [F_y, mu_y] = mFormulaLateral(alpha, gamma, F_z, p)
%Pacejka lateral tire forces (pure slip)

%turn slip is neglected, so xi_i=1;
%all scaling factors lambda = 1;

%coordinate system transformation
%alpha = -alpha;

S_hy = sign(gamma)*(p.p_hy1 + p.p_hy3*abs(gamma));
S_vy = sign(gamma)*F_z*(p.p_vy1 + p.p_vy3*abs(gamma));

alpha_y = alpha + S_hy;
mu_y = p.p_dy1*(1-p.p_dy3*gamma^2);

C_y = p.p_cy1;
D_y = mu_y*F_z;
E_y = p.p_ey1;
K_y = F_z*p.p_ky1; %simplify K_y0 to p.p_ky1*F_z
B_y = K_y/(C_y*D_y);

%magic tire formula
F_y = D_y*sin(C_y*atan(B_y*alpha_y - E_y*(B_y*alpha_y - atan(B_y*alpha_y)))) + S_vy;

end

function F_y = mFormulaLateralComb(kappa, alpha, gamma, mu_y, F_z, F0_y, p)
%lateral tire forces (combined slip)

%turn slip is neglected, so xi_i=1;
%all scaling factors lambda = 1;

S_hykappa = p.r_hy1; 

kappa_s = kappa + S_hykappa;

B_ykappa = p.r_by1*cos(atan(p.r_by2*(alpha-p.r_by3)));
C_ykappa = p.r_cy1;
E_ykappa = p.r_ey1;
D_ykappa = F0_y/(cos(C_ykappa*atan(B_ykappa*S_hykappa - E_ykappa*(B_ykappa*S_hykappa - atan(B_ykappa*S_hykappa)))));

D_vykappa = mu_y*F_z*(p.r_vy1 + p.r_vy3*gamma) * cos(atan(p.r_vy4*alpha));
S_vykappa = D_vykappa*sin(p.r_vy5*atan(p.r_vy6*kappa));

%magic tire formula
F_y = D_ykappa*cos(C_ykappa*atan(B_ykappa*kappa_s - E_ykappa*(B_ykappa*kappa_s - atan(B_ykappa*kappa_s)))) + S_vykappa;
end

function F_x = mFormulaLongitudinal(kappa, gamma, F_z, p)
%Pacejka longitudinal tire forces (pure slip)

%turn slip is neglected, so xi_i=1;
%all scaling factors lambda = 1;

%coordinate system transformation
kappa = -kappa;

S_hx = p.p_hx1;
S_vx = F_z*p.p_vx1;

kappa_x = kappa + S_hx;
mu_x = p.p_dx1*(1-p.p_dx3*gamma^2);

C_x = p.p_cx1;
D_x = mu_x*F_z;
E_x = p.p_ex1;
K_x = F_z*p.p_kx1;
B_x = K_x/(C_x*D_x);

%magic tire formula
F_x = D_x*sin(C_x*atan(B_x*kappa_x - E_x*(B_x*kappa_x - atan(B_x*kappa_x))) + S_vx);
end

function F_x = mFormulaLongitudinalComb(kappa, alpha, F0_x, p)
%Pacejka longitudinal tire forces (combined slip)

%turn slip is neglected, so xi_i=1;
%all scaling factors lambda = 1;

S_hxalpha = p.r_hx1; 

alpha_s = alpha + S_hxalpha;

B_xalpha = p.r_bx1*cos(atan(p.r_bx2*kappa));
C_xalpha = p.r_cx1;
E_xalpha = p.r_ex1;
D_xalpha = F0_x/(cos(C_xalpha*atan(B_xalpha*S_hxalpha - E_xalpha*(B_xalpha*S_hxalpha - atan(B_xalpha*S_hxalpha)))));

%magic tire formula
F_x = D_xalpha*cos(C_xalpha*atan(B_xalpha*alpha_s - E_xalpha*(B_xalpha*alpha_s - atan(B_xalpha*alpha_s))));
end



%------------- END OF CODE --------------
