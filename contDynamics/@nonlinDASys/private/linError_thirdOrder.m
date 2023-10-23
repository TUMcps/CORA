function [error, errorInt, errorInt_x, errorInt_y, R_y] = linError_thirdOrder(obj, options, R, Verror_y)
% linError_thirdOrder - computes the linearization error using a third
% order Taylor expansion
%
% Syntax:
%    [error, errorInt, errorInt_x, errorInt_y, R_y] = ...
%           linError_thirdOrder(obj, options, R, Verror_y)
%
% Inputs:
%    obj - nonlinear differential algebraic system object
%    options - options struct
%    R - actual reachable set
%    Verror_y - set of algebraic linearization error
%
% Outputs:
%    error - zonotope overapproximating the linearization error
%    errorInt - interval overapproximating the linearization error
%    errorInt_x - interval overapproximating the linearization error (dynamic part)
%    errorInt_y - interval overapproximating the linearization error (constraint part)
%    R_y - reachable set of the algebraic part
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       21-June-2013
% Last update:   16-June-2016
%                25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set handle to correct file
obj = setHessian(obj,'standard');
obj = setThirdOrderTensor(obj,'int');

%compute set of algebraic variables
f0_con = obj.linError.f0_con;
D = obj.linError.D;
E = obj.linError.E;
F_inv = obj.linError.F_inv;
R_y_cor = -F_inv*(f0_con + D*R); %correlated part
R_y_add = -F_inv*(E*options.U + Verror_y); %uncorrelated part


%obtain intervals and combined interval z
dx = interval(R);
dy = interval(R_y_cor + R_y_add);
du = interval(options.U);
dz = [dx; dy; du];

%compute interval of reachable set
totalInt_x = dx + obj.linError.p.x;

%compute interval of algebraic states
totalInt_y = dy + obj.linError.p.y;

%compute intervals of input
totalInt_u = du + obj.linError.p.u;

%obtain hessian and third order tensor
[Hf, Hg] = obj.hessian(obj.linError.p.x, obj.linError.p.y, obj.linError.p.u);
[Tf, Tg] = obj.thirdOrderTensor(totalInt_x, totalInt_y, totalInt_u);

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
R_xyu = cartProd(R_xy, options.U);
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
Z_err_x_add = get(obj.linError.CF_inv*error_y_secondOrder,'Z');
error_secondOrder = zonotope(Z_err_x + Z_err_x_add);
error_thirdOrder = error_thirdOrder_x_zono + obj.linError.CF_inv*error_thirdOrder_y_zono;
error = error_secondOrder + error_thirdOrder;

%reduce
error = reduce(error,options.reductionTechnique,options.zonotopeOrder);
error_y = reduce(error_y,options.reductionTechnique,options.zonotopeOrder);

%update R_y
R_y =  obj.linError.p.y + (-F_inv)*(f0_con + D*R + E*options.U + error_y);

%error intervals
errorIHabs = abs(interval(error));
errorInt = supremum(errorIHabs);

errorIHabs_y = abs(interval(error_y));
errorInt_y = supremum(errorIHabs_y);

errorIHabs_x = abs(interval(error_x));
errorInt_x = supremum(errorIHabs_x);

% ------------------------------ END OF CODE ------------------------------
