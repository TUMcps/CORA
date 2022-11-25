function [obj,A_lin,U] = linearize(obj,R,options)
% linearize - linearizes the nonlinearSysDT object
%
% Syntax:  
%    [obj,A_lin,U] = linearize(obj,R,options)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    R - initial reachable set
%    options - options struct
%
% Outputs:
%    obj - nonlinearSysDT system object with additional properties
%    A_lin - system matrix of the linearized system
%    U - reachable set due to the inputs
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  29-January-2018
% Last revision:---

%------------- BEGIN CODE --------------

%linearization point p.u of the input is the center of the input u
p.u = center(options.U) + options.uTrans;

%linearization point p.x and p.y
x0 = center(R);
p.x = x0;

%substitute p into the system equation in order to obtain the constant
%input
f0 = obj.mFile(p.x, p.u);

%get jacobian matrices
[A_lin,B_lin] = obj.jacobian(p.x, p.u);


uTrans = f0; %B*Ucenter from linOptions.U not added as the system is linearized around center(U)
Udelta = B_lin*(options.U+(-center(options.U)));
U = Udelta + uTrans;

%save linearization point
obj.linError.p=p;

%------------- END OF CODE --------------