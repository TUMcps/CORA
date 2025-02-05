function [nlnsysDA,linsys,linParams,linOptions] = priv_linearize(nlnsysDA,R,R_y,params,options)
% priv_linearize - linearizes the nonlinear differential algebraic system;
%    linearization error is not included (see linError)
%
% Syntax:
%    [nlnsysDA,linsys,linParams,linOptions] = priv_linearize(nlnsysDA,R,R_y,params,options)
%
% Inputs:
%    nlnsysDA - nonlinear DAE system object
%    R - 
%    R_y - 
%    params - model parameters
%    options - options struct
%
% Outputs:
%    nlnsysDA - nonlinDASys object
%    linsys - linear system object
%    linParams - model parameters for the linearized system
%    linOptions - options for the linearized system
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       21-November-2011
% Last update:   23-May-2013       
%                28-July-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%linearization point p.u of the input is the center of the input u
p.u = center(params.U) + params.uTrans;

%linearization point p.x and p.y
x0 = center(R);
y0 = center(R_y);
f0prev_x = nlnsysDA.dynFile(x0, y0, p.u);


try %if time step already created
    p.x = x0 + f0prev_x*0.5*options.timeStep;
    p.y = aux_consistentState(nlnsysDA, p.x, y0, p.u);
catch
    disp('time step not yet created; this message should only appear once!');
    p.x = x0;
    p.y = y0;
end

%substitute p into the system equation in order to obtain the constant
%input
f0_dyn = nlnsysDA.dynFile(p.x, p.y, p.u);
f0_con = nlnsysDA.conFile(p.x, p.y, p.u);

%get jacobian matrices
[A,B,C,D,E,F] = nlnsysDA.jacobian(p.x, p.y, p.u);


%compute matrices of the linearized system
F_inv = pinv(F);
CF_inv = C*F_inv;
f0 = f0_dyn - CF_inv*f0_con;
A_lin = A - CF_inv*D;
B_lin = B - CF_inv*E;
max(real(eig(A_lin)));


% set up params and options for linearized system
linOptions=options;
linParams.uTrans = f0; %B*Ucenter from linOptions.U not added as the system is linearized around center(U)
linParams.U = B_lin*(params.U+(-center(params.U)));
linOptions.originContained = false;

%set up linearized system
linsys = linearSys('linsys',A_lin,1); %B=1 as input matrix encountered in uncertain inputs

%save constant input and matrices
nlnsysDA.linError.f0 = f0;
nlnsysDA.linError.f0_con = f0_con;
nlnsysDA.linError.D = D;
nlnsysDA.linError.E = E;
nlnsysDA.linError.F_inv = F_inv;
nlnsysDA.linError.CF_inv = CF_inv;

%save linearization point
nlnsysDA.linError.p=p;

end


% Auxiliary functions -----------------------------------------------------

function y0 = aux_consistentState(nlnsysDA, x0, y0, u0)

while true
    l = nlnsysDA.conFile(x0, y0, u0);
    [~,~,~,~,~,F] = nlnsysDA.jacobian(x0, y0, u0);

    %evaluate jacobian
    delta_y = F\(-l);
    
    %check convergence
    if norm(delta_y)<1e-10
        break
    end
    
    %update steady state solution
    y0 = y0 + delta_y;
end

end

% ------------------------------ END OF CODE ------------------------------
