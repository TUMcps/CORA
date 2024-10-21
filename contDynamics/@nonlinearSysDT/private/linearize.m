function [nlnsysDT,A_lin,U] = linearize(nlnsysDT,R,params)
% linearize - linearizes the nonlinearSysDT object
%
% Syntax:
%    [nlnsysDT,A_lin,U] = linearize(nlnsysDT,R,params)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    R - initial reachable set
%    params - model parameters
%
% Outputs:
%    obj - nonlinearSysDT system object with additional properties
%    A_lin - system matrix of the linearized system
%    U - reachable set due to the inputs
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       21-August-2012
% Last update:   29-January-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%linearization point p.u of the input is the center of the input set U
p.u = center(params.U);

%linearization point p.x and p.y
x0 = center(R);
p.x = x0;

%substitute p into the system equation in order to obtain the constant
%input
f0 = nlnsysDT.mFile(p.x, p.u);

%get jacobian matrices
[A_lin,B_lin] = nlnsysDT.jacobian(p.x, p.u);

Udelta = B_lin*(params.U+(-center(params.U))); %B*Ucenter from linOptions.U not added as the system is linearized around center(U)
U = Udelta + f0;

%save linearization point
nlnsysDT.linError.p=p;

% ------------------------------ END OF CODE ------------------------------
