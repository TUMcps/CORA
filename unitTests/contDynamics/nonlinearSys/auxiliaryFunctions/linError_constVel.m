function [error] = linError_constVel(obj,w,IH,IH_y)
% linError - computes the linearization error
%
% Syntax:  
%    [obj] = linError(obj,options)
%
% Inputs:
%    obj - nonlinear system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    obj - nonlinear system object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      15-March-2012 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% %compute intervals of total reachable set
% totalInt=interval(IH);

%obtain model factors
c = modelFactors(w);
c_abs = abs(c);

%helper values
min_x = infimum(IH);
max_x = supremum(IH);

%extreme absolue values 
%x_1
max_abs_x1 = max(abs(min_x(1)),abs(max_x(1))); 
%x_3
max_abs_x3 = max(abs(min_x(3)),abs(max_x(3))); 
%x_4
min_abs_x4 = 7.4; %hard coded!!
max_abs_x4 = 7.6; %hard coded!!

%auxiliary variables
aux_min_1 = w(1) - infimum(IH(4)) - infimum(IH_y(1));
aux_max_1 = w(1) - supremum(IH(4)) - supremum(IH_y(1));
tilde_IM_1 = max(abs(aux_min_1),abs(aux_max_1));

aux_min_2 = w(2) - infimum(IH(5)) - infimum(IH_y(2));
aux_max_2 = w(2) - supremum(IH(5)) - supremum(IH_y(2));
tilde_IM_2 = max(abs(aux_min_2),abs(aux_max_2));

aux_min_3 = w(3) - infimum(IH(2)) - infimum(IH_y(3));
aux_max_3 = w(3) - supremum(IH(2)) - supremum(IH_y(3));
tilde_IM_3 = max(abs(aux_min_3),abs(aux_max_3));

aux_min_4 = w(4) - infimum(IH(3)) - infimum(IH_y(4));
aux_max_4 = w(4) - supremum(IH(3)) - supremum(IH_y(4));
tilde_IM_4 = max(abs(aux_min_4),abs(aux_max_4));


%cos and sin max
min_angle = min_x(1) + min_x(2);
max_angle = max_x(1) + max_x(2);

%ASSUME that angles are between -pi and pi!!
%cosinus ranges
if max_angle<0
    cos_min = cos(min_angle);
    cos_max = cos(max_angle);
elseif min_angle>0
    cos_min = cos(max_angle);
    cos_max = cos(min_angle);
else
    cos_min = min(cos(min_angle), cos(max_angle));
    cos_max = 1;
end

%sinus ranges
if max_angle<-pi/2
    sin_min = sin(max_angle);
    sin_max = sin(min_angle);
elseif min_angle>pi/2
    sin_min = sin(max_angle);
    sin_max = sin(min_angle);
elseif min_angle>-pi/2 && max_angle<pi/2
    sin_min = sin(min_angle);
    sin_max = sin(max_angle);
elseif min_angle<-pi/2 && max_angle<pi/2
    sin_min = -1;
    sin_max = max(sin(min_angle),sin(max_angle));
elseif min_angle>-pi/2 && max_angle>pi/2
    sin_min = min(sin(min_angle),sin(max_angle));
    sin_max = 1;
end

%auxiliary variables for sin and cos
Sm = max(abs(sin_min),abs(sin_max));
Cm = max(abs(cos_min),abs(cos_max));

%translate intervals by linearization point
IH=IH+(-obj.linError.p.x);

%obtain maximum absolute values within IH, IHinput
IHinf=abs(infimum(IH));
IHsup=abs(supremum(IH));
dx_new=max(IHinf,IHsup);

%conversion
dx(4,1) = 0.2; %delta v; hard coded!!
dx(5,1) = dx_new(4,1);
dx(6,1) = dx_new(5,1);

IHinf_y=abs(infimum(IH_y));
IHsup_y=abs(supremum(IH_y));
dy=max(IHinf_y,IHsup_y);


%compute linearization error by hand-derived Lagrange remainders
%first coordinate
error(1,1) = 2/(min_abs_x4^2)*(c_abs(2)*dx(1) + c_abs(5)*(dx(2) + dy(3)) ...
    + c_abs(3)*(dx(5) + dy(1)) + c_abs(4)*(dx(6) + dy(2)))*dx(4) ...
    + 2/(min_abs_x4^3)*(c_abs(2)*max_abs_x1 + c_abs(3)*tilde_IM_1 ...
    + c_abs(4)*tilde_IM_2 + c_abs(5)*tilde_IM_3 + c_abs(6)*tilde_IM_4)*dx(4)^2 ...
    + 4/(min_abs_x4^3)*c_abs(1)*dx(3)*dx(4) ...
    + 6/(min_abs_x4^4)*c_abs(1)*max_abs_x3*dx(4)^2;

error(2,1) = 0;

error(3,1) = 2/(min_abs_x4^2)*c_abs(7)*dx(3)*dx(4) ...
    + 2/(min_abs_x4^3)*c_abs(7)*max_abs_x3*dx(4)^2;

error(4,1) = max_abs_x4*Cm*(dx(1) + dx(2))^2 + 2*Sm*(dx(1) + dx(2))*dx(4);
error(5,1) = max_abs_x4*Sm*(dx(1) + dx(2))^2 + 2*Cm*(dx(1) + dx(2))*dx(4);

error(6,1) = 0;

%divide by 2 since the second order Taylor term has factor 1/2!
error = 0.5*error;

%------------- END OF CODE --------------
