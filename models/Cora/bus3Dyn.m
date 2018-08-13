function f=bus3Dyn(t,x,y,u)

%x(1) = omega
%x(2) = T_m

%y(1) = E
%y(2) = V_2
%y(3) = V_3
%y(4) = Theta_1 - delta
%y(5) = Theta_2 - delta
%y(6) = Theta_3 - delta

%u(1) = P_c
%u(2) = P_w

%slack bus voltage
V_1 = 1;

%parameters
X_m = 0.2;
M = 1/(15*pi);
D = 0.04;
T_sv = 1; %this was wrongly referred to as T_m in the table
omega_s = 120*pi;
R_d = 0.05; %guessed value!


%dynamic model
f(1,1) = -D/M*x(1) + 1/M*x(2) - (y(1)*V_1)/(M*X_m)*sin(-y(4)) + D/M*omega_s; %omega
f(2,1) = -1/(T_sv*R_d*omega_s)*x(1) - 1/T_sv*x(2) + 1/(T_sv*R_d) + 1/T_sv*u(1); %T_m

