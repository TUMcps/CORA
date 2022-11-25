function [E_err, c_tilde, d_tilde] = trajectoryError(obj,f,u,R0,tmin,tmax,tc,order,n)
% trajectoryError - computes the error between a linear system solution and
% a constant derivative solution. This function is required for computing
% the error for the mapping-based approach for guard intersection as
% presented at HSCC'12.
%
% Syntax:  
%    E_err = trajectoryError(A,f,u,R0,tmin,tmax,tc,order,n)
%
% Inputs:
%    obj - linear system object
%    f - constant flow offset of other dynamics
%    u - input of linear system
%    R0 - initial set
%    tmin - start time
%    tmax - final time
%    tc - "center time", i.e. hitting time for x0
%    
%
% Outputs:
%    E_err - abstraction error
%
% Example: 
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      23-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%in case there is no flow correction
if ~iscell(f)
    fTmp = f;
    f = [];
    f{1} = fTmp;
    f{2} = fTmp;
end

%obtain dimension and system matrix
dim = length(u);
A = obj.A;

%compute powers of tmin and tmax
tZono = timePowers(tmin,tmax,order);

%obtain center and delta set
x0 = center(R0);
x0zono = zonotope(x0);
uzono = zonotope(u);
Y = R0 + (-x0);

%auxiliary vector
E_err = 0;
E_err_u = 0;
u_tilde = 0;
for i=2:order
    E_err = E_err + A^i/factorial(i)*(tZono{i-1}*Y + (tZono{i-1} + (-1)*tc^(i-1))*x0zono);
    u_tilde = u_tilde + A^i/factorial(i)*(tZono{1} * (-tc^(i-1)))*x0zono;
end
for i=1:order
    E_err_u = E_err_u + A^i/factorial(i+1)*(tZono{i} + (-1)*tc^(i))*uzono;
end
E_err = E_err + E_err_u;
E_err = tZono{1}*E_err;
u_tilde = u_tilde + tZono{1}*E_err_u;

%compute remainder matrix
A_abs = abs(A);
Apower_abs = eye(length(A_abs));
M = eye(length(A_abs));

for i=1:order
    Apower_abs = Apower_abs*A_abs;
    M = M + Apower_abs*tmax^i/factorial(i);
end 
  
%determine error due to finite Taylor series
W = expm(A_abs*tmax) - M;
W = abs(W); % for computational reasons when numbers are very small
%compute absolute value of W for numerical stability
E = interval(-W,W);

%consider E
E_err = E_err + E*R0 + E*tmax*uzono + tZono{1}*zonotope(f{1}-f{2}); %error due to flow correction
u_tilde = u_tilde + E*tmax*uzono + tZono{1}*zonotope(f{1}-f{2});

%obtain c_tilde and d_tilde
A_tilde = zeros(dim);
for i=2:order
    A_tilde = A_tilde + A^i/factorial(i)*tZono{i};
end
A_tilde = A_tilde + E;
c_tilde = A_tilde.'*n;
c_tilde = intervalMatrix(c_tilde);
c_tilde = c_tilde.int;

d_tilde = -n.'*u_tilde;
d_tilde = interval(d_tilde);




function tZono = timePowers(tmin,tmax,order)


%first order
tminPow(1) = tmin;
tmaxPow(1) = tmax;

tradPow(1) = 0.5*(tmaxPow(1) - tminPow(1));
tmidPow(1) = 0.5*(tmaxPow(1) + tminPow(1));

tRad{1} = tradPow(1);
tZono{1} = matZonotope(tmidPow(1),tRad);

%higher orders
for i=1:order
    tminPow(i+1) = tminPow(i)*tmin;
    tmaxPow(i+1) = tmaxPow(i)*tmax;
    
    tradPow(i+1) = 0.5*(tmaxPow(i+1) - tminPow(i+1));
    tmidPow(i+1) = 0.5*(tmaxPow(i+1) + tminPow(i+1));
    
    tRad{1} = tradPow(i+1);
    tZono{i+1} = matZonotope(tmidPow(i+1),tRad);
end




%------------- END OF CODE --------------