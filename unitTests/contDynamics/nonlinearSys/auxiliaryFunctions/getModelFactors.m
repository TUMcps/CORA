function getModelFactors(p,mu_max)
% getModelFactors - obtains model factors for computing the linearization
% error
%
% Syntax:  
%    getModelFactors(p,mu_max)
%
% Inputs:
%    p - vehicle parameters
%    mu_max - maximum tire-road friction
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
% Written:      04-May-2011
% Last update:  31-August-2011
% Last revision:---

%------------- BEGIN CODE --------------


%create equivalent bicycle parameters
g = 9.81; %[m/s^2]
mu = mu_max;
C_Sf = -p.tire.p_ky1/p.tire.p_dy1; 
C_Sr = -p.tire.p_ky1/p.tire.p_dy1; 
lf = p.a;
lr = p.b;
%h = p.h_s;
h = 0;
m = p.m;
I = p.I_z;

syms uL3R uL6R

%obtain control gains
[carInput,k] = DOTcontrol(zeros(6,1),zeros(5,1),zeros(5,1));

%model factors
c(1,1) =  mu/(lr+lf)*(C_Sr*(g*lf + uL6R*h)*lr - C_Sf*(g*lr - uL6R*h)*lf); %c_doc_1
c(2,1) = -mu/(lr+lf)*(C_Sr*(g*lf + uL6R*h) + C_Sf*(g*lr-uL6R*h)); %c_doc_2
c(3,1) = -mu/(lr+lf)*(C_Sf*(g*lr-uL6R*h))*k(1)*sin(uL3R); %c_doc_3
c(4,1) =  mu/(lr+lf)*(C_Sf*(g*lr-uL6R*h))*k(1)*cos(uL3R); %c_doc_4
c(5,1) =  mu/(lr+lf)*(C_Sf*(g*lr-uL6R*h))*k(2); %c_doc_5
c(6,1) =  mu/(lr+lf)*(C_Sf*(g*lr-uL6R*h))*k(3); %c_doc_6
c(7,1) = -mu*m/(I*(lr+lf))*(lf^2*C_Sf*(g*lr-uL6R*h) + lr^2*C_Sr*(g*lf + uL6R*h)); %c_doc_7

%write mfile for model factors
createModelFactors(c);

%------------- END OF CODE --------------