function [obj,linSys,options,linOptions] = linearize(obj,options,R,R_y)
% linearize - linearizes the nonlinear differential algebraic system;
% linearization error is not included (see linError)
%
% Syntax:  
%    [obj,linSys,options,linOptions] = linearize(obj,options,R)
%
% Inputs:
%    obj - nonlinear DAE system object
%    options - options struct
%
% Outputs:
%    obj - linear system object
%    linSys - linear system object
%    options - options struct
%    linOptions - options for the linearized system
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      21-November-2011
% Last update:  23-May-2013       
%               28-July-2020
% Last revision:---

%------------- BEGIN CODE --------------

%linearization point p.u of the input is the center of the input u
p.u = center(options.U) + options.uTrans;

%linearization point p.x and p.y
x0 = center(R);
y0 = center(R_y);
f0prev_x = obj.dynFile(x0, y0, p.u);


try %if time step already created
    p.x = x0 + f0prev_x*0.5*options.timeStep;
    p.y = consistentState(obj, p.x, y0, p.u);
catch
    disp('time step not yet created; this message should only appear once!');
    p.x = x0;
    p.y = y0;
end

%substitute p into the system equation in order to obtain the constant
%input
f0_dyn = obj.dynFile(p.x, p.y, p.u);
f0_con = obj.conFile(p.x, p.y, p.u);

%get jacobian matrices
[A,B,C,D,E,F] = obj.jacobian(p.x, p.y, p.u);


%compute matrices of the linearized system
F_inv = pinv(F);
CF_inv = C*F_inv;
f0 = f0_dyn - CF_inv*f0_con;
A_lin = A - CF_inv*D;
B_lin = B - CF_inv*E;
max(real(eig(A_lin)));


%set up otions for linearized system
linOptions=options;

linOptions.uTrans = f0; %B*Ucenter from linOptions.U not added as the system is linearized around center(U)
linOptions.U = B_lin*(options.U+(-center(options.U)));
linOptions.originContained = 0;

%set up linearized system
linSys = linearSys('linSys',A_lin,1); %B=1 as input matrix encountered in uncertain inputs

%save constant input and matrices
obj.linError.f0 = f0;
obj.linError.f0_con = f0_con;
obj.linError.D = D;
obj.linError.E = E;
obj.linError.F_inv = F_inv;
obj.linError.CF_inv = CF_inv;

%save linearization point
obj.linError.p=p;


function y0 = consistentState(obj, x0, y0, u0)

%init
converged = 0;

while ~converged
    
    l = obj.conFile(x0, y0, u0);
    [~,~,~,~,~,F] = obj.jacobian(x0, y0, u0);


    %evaluate jacobian
    delta_y = F\(-l);
    
    %check convergence
    if norm(delta_y)<1e-10
        converged = 1;
    end
    
    %update steady state solution
    y0 = y0 + delta_y;
end


%------------- END OF CODE --------------