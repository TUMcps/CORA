function [error,errorInt,errorInt_x,errorInt_y,R_y] = linError_thirdOrder(nlnsysDA,R,Verror_y,params,options)
% linError_thirdOrder - computes the linearization error using a third
%    order Taylor expansion
%
% Syntax:
%    [error,errorInt,errorInt_x,errorInt_y,R_y] = ...
%           linError_thirdOrder(obj,R,Verror_y,params,options)
%
% Inputs:
%    nlnsysDA - nonlinear differential algebraic system object
%    R - actual reachable set
%    Verror_y - set of algebraic linearization error
%    params - model parameters
%    options - options struct
%
% Outputs:
%    error - zonotope overapproximating the linearization error
%    errorInt - interval overapproximating the linearization error
%    errorInt_x - interval overapproximating the linearization error (dynamic part)
%    errorInt_y - interval overapproximating the linearization error (constraint part)
%    R_y - reachable set of the algebraic part
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       21-June-2013
% Last update:   16-June-2016
%                25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set handle to correct file
nlnsysDA = setHessian(nlnsysDA,'standard');
nlnsysDA = setThirdOrderTensor(nlnsysDA,'int');

%compute set of algebraic variables
f0_con = nlnsysDA.linError.f0_con;
D = nlnsysDA.linError.D;
E = nlnsysDA.linError.E;
F_inv = nlnsysDA.linError.F_inv;
R_y_cor = -F_inv*(f0_con + D*R); %correlated part
R_y_add = -F_inv*(E*params.U + Verror_y); %uncorrelated part


%obtain intervals and combined interval z
dx = interval(R);
dy = interval(R_y_cor + R_y_add);
du = interval(params.U);
dz = [dx; dy; du];

%compute interval of reachable set
totalInt_x = dx + nlnsysDA.linError.p.x;

%compute interval of algebraic states
totalInt_y = dy + nlnsysDA.linError.p.y;

%compute intervals of input
totalInt_u = du + nlnsysDA.linError.p.u;

%obtain hessian and third order tensor
[Hf, Hg] = nlnsysDA.hessian(nlnsysDA.linError.p.x, nlnsysDA.linError.p.y, nlnsysDA.linError.p.u);
[Tf, Tg] = nlnsysDA.thirdOrderTensor(totalInt_x, totalInt_y, totalInt_u);

%store Hf and Hg as real-valued 
for i=1:length(Hf)
    Hf{i} = center(Hf{i});
end
for i=1:length(Hg)
    Hg{i} = center(Hg{i});
end


%compute zonotope of state, constarint variables, and input
Z_x = [R.c,R.G];
Z_y_cor = [R_y_cor.c,R_y_cor.G];
Z_y_add = [R_y_add.c,R_y_add.G];
Z_0 = zeros(length(Z_x(:,1)), length(Z_y_add(1,:)));
R_xy = zonotope([Z_x, Z_0; Z_y_cor, Z_y_add]);
R_xyu = cartProd(R_xy, params.U);
R_xyu = reduce(R_xyu,options.reductionTechnique,options.errorOrder);


%obtain absolute values
dz_abs = max(abs(infimum(dz)), abs(supremum(dz)));

%second order
error_x_secondOrder = 0.5*quadMap_parallel(R_xyu, Hf);
error_y_secondOrder = 0.5*quadMap_parallel(R_xyu, Hg);

%third order interval evaluation (dynamic part)
for i=1:length(Tf(:,1))
    error_sum = interval(0,0);
    for j=1:length(Tf(1,:))
        error_tmp = dz'*Tf{i,j}*dz;
        error_sum = error_sum + error_tmp * dz(j);
    end
    error_x_thirdOrder(i,1) = 1/6*error_sum;
end

%third order interval evaluation (algebraic part)
for i=1:length(Tg(:,1))
    error_sum = interval(0,0);
    for j=1:length(Tg(1,:))
        error_tmp = dz'*Tg{i,j}*dz;
        error_sum = error_sum + error_tmp * dz(j);
    end
    error_y_thirdOrder(i,1) = 1/6*error_sum;
end

%convert to zonotopes
error_thirdOrder_x_zono = zonotope(error_x_thirdOrder);
error_thirdOrder_y_zono = zonotope(error_y_thirdOrder);

%combine results
error_x = error_x_secondOrder + error_thirdOrder_x_zono;
error_y = error_y_secondOrder + error_thirdOrder_y_zono;

%compute final error
Z_err_x = [error_x_secondOrder.c,error_x_secondOrder.G];
Z_err_x_add = nlnsysDA.linError.CF_inv*error_y_secondOrder;
Z_err_x_add_mat = [Z_err_x_add.c, Z_err_x_add.G];
error_secondOrder = zonotope(Z_err_x + Z_err_x_add_mat);
error_thirdOrder = error_thirdOrder_x_zono + nlnsysDA.linError.CF_inv*error_thirdOrder_y_zono;
error = error_secondOrder + error_thirdOrder;

%reduce
error = reduce(error,options.reductionTechnique,options.zonotopeOrder);
error_y = reduce(error_y,options.reductionTechnique,options.zonotopeOrder);

%update R_y
R_y =  nlnsysDA.linError.p.y + (-F_inv)*(f0_con + D*R + E*params.U + error_y);

%error intervals
errorIHabs = abs(interval(error));
errorInt = supremum(errorIHabs);

errorIHabs_y = abs(interval(error_y));
errorInt_y = supremum(errorIHabs_y);

errorIHabs_x = abs(interval(error_x));
errorInt_x = supremum(errorIHabs_x);

% ------------------------------ END OF CODE ------------------------------
