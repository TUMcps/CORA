function res = test_linear_trajectoryError()
% test_linear_trajectoryError - unit_test_function of trajectoryError
% function required for the mapping-based guard intersection as described
% in HSCC'12.
%
% Syntax:  
%    res = test_linear_trajectoryError()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      23-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% obtain linear system
A = [-1 -4; 4 -1];
linSys = linearSys('abstractionErrorSys',A,1);

% constant flow offset of other dynamics
f = [0; -0.3];

% input of linear system
u = [0; 1];

% initial set
R0 = zonotope([0.6 0.1 0; 1 0 0.1]);

%  start time
tmin = 0;

% final time
tmax = 0.05;

% "center time", i.e. hitting time for x0
tc = 0.025;

% Taylor order
order = 4;

% halfspace normal
n = [1; 0];

% call function
E_err = trajectoryError(linSys,f,u,R0,tmin,tmax,tc,order,n);

% obtain values for comparison---------------------------------------------
% powers of time intervals
tInt = interval(tmin,tmax);

for i = 1:(order+1)
    res = tInt^i;
    c = mid(res);
    gen{1} = rad(res);
    tZono{i} = matZonotope(c,gen); 
end

%separation of initial uncertainty
xi = zonotope(center(R0));
Y = R0 + (-1*xi);

% state error
E_state = 0;
for i=2:order
    E_state = E_state + A^i/factorial(i)*(tZono{i-1}*Y + (tZono{i-1} + (-1)*tc^(i-1))*xi);
end

% input error
E_input = 0;
for i=1:order
    E_input = E_input + A^i/factorial(i+1)*(tZono{i} + (-1)*tc^(i))*zonotope(u);
end

% state remainder error
Aabs = abs(A);
Waux = 0;
for i=0:order
    Waux = Waux + Aabs^i*tmax^i/factorial(i);
end
W = expm(Aabs*tmax) - Waux;
Ehat = interval(-W,W);

% total error
E_err_calc = tZono{1}*(E_state + E_input) + Ehat*R0 + Ehat*tmax*zonotope(u);

% enclose results by intervals
I_err = interval(E_err);
I_err_calc = interval(E_err_calc);

% check if slightly bloated versions enclose each other
res_1 = (I_err <= enlarge(I_err_calc,1+1e-8));
res_2 = (I_err_calc <= enlarge(I_err,1+1e-8));

% final result
res = res_1*res_2;

%--------------------------------------------------------------------------

%------------- END OF CODE --------------
