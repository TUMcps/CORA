function nlnsysDT = nonlinearSysDT(nlnARX)
% nonlinearSysDT - convert nonlinearARX sys to nonlinearSysDT
%
% Syntax:
%    nlnsysDT = nonlinearSysDT(nlnARX)
%
% Inputs:
%    nlnARX - nonlinearARX object
%
% Outputs:
%    nlnsysDT - nonlinearSysDT object
%
% Example:
%     f = @(y,u) [0.5*y(1,1) + u(1,1) - cos(u(2,1)); ...
%         0.4*y(3,1) + u(2,1)*cos(y(1,1)); 0.6*y(5,1) + u(4,1)*sin(y(1,1))];
%     dt = 0.1; dim_y = 3; dim_u = 2; n_p = 2;
%     nlnARX = nonlinearARX(f,dt,dim_y,dim_u,n_p);
%     nlnsysDT = nonlinearSysDT(nlnARX);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT, nonlinearARX

% Authors:       Laura Luetzow
% Written:       13-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

states = nlnARX.nrOfOutputs * nlnARX.n_p;
inputs = nlnARX.nrOfInputs * (nlnARX.n_p + 1);
name = sprintf('%s_conv', nlnARX.name);
fun = @(x,u) [x(nlnARX.nrOfOutputs+1:end); nlnARX.mFile(x,u)];
nlnsysDT = nonlinearSysDT(name,fun,nlnARX.dt,states,inputs);

% ------------------------------ END OF CODE ------------------------------
