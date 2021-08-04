function [error, errorInt, errorInt_x, errorInt_y, R_y] = linError_mixed_noInt(obj, options, R, Verror_y)
% linError_mixed_noInt - computes the linearization error
%
% Syntax:  
%    [error, errorInt, errorInt_x, errorInt_y, R_y] = ...
%           linError_mixed_noInt(obj, options, R, Verror_y)
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

% Author:       Matthias Althoff
% Written:      21-November-2011
% Last update:  23-May-2013
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

% set handle to correct file
obj = setHessian(obj,'standard');

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

%obtain hessian tensor
[Hf, Hg] = obj.hessian(totalInt_x, totalInt_y, totalInt_u);

% %[Hf_old, Hg_old] = hessianTensor(totalInt_x, totalInt_y, totalInt_u);
% %[Hf, Hg] = hessianTensor_optimized(totalInt_x, totalInt_y, totalInt_u);
% if strcmp(options.mode,'normal')
%     [Hf, Hg] = hessianTensor_normal_monotone_parallel(totalInt_x, totalInt_y, totalInt_u);
%     %[Hf, Hg] = hessianTensor_normal_parallel(totalInt_x, totalInt_y, totalInt_u);
%     %[Hf,Hg]=hessianTensor_normal_test2;
%     %[Hf, Hg] = hessianTensor_normal_test(totalInt_x, totalInt_y, totalInt_u);
%     %[Hf, Hg] = hessianTensor_normal(totalInt_x, totalInt_y, totalInt_u);
% elseif strcmp(options.mode,'fault')
%     [Hf, Hg] = hessianTensor_fault_monotone_parallel(totalInt_x, totalInt_y, totalInt_u);
%     %[Hf, Hg] = hessianTensor_fault(totalInt_x, totalInt_y, totalInt_u);
% end

%compute zonotope of state, constarint variables, and input
Z_x = R.Z;
Z_y_cor = R_y_cor.Z;
Z_y_add = R_y_add.Z;
Z_0 = zeros(length(Z_x(:,1)), length(Z_y_add(1,:)));
R_xy = zonotope([Z_x, Z_0; Z_y_cor, Z_y_add]);
R_xyu = cartProd(R_xy, options.U);
R_xyu = reduce(R_xyu,options.reductionTechnique,options.errorOrder);

%obtain absolute values
dz_abs = max(abs(infimum(dz)), abs(supremum(dz)));
%dz_abs = max(abs(inf(dz)), abs(sup(dz))); %INTLAB syntax

%separate evaluation
for i=1:length(Hf)
    Hf_mid{i} = sparse(center(Hf{i}));
    Hf_rad{i} = sparse(rad(Hf{i}));
end
for i=1:length(Hg)
    Hg_mid{i} = sparse(center(Hg{i}));
    Hg_rad{i} = sparse(rad(Hg{i}));
end
%zonotope evaluation
% error_x_mid_old = 0.5*quadMap(R_xyu, Hf_mid);
% error_y_mid_old = 0.5*quadMap(R_xyu, Hg_mid);
error_x_mid = 0.5*quadMap_parallel(R_xyu, Hf_mid);
error_y_mid = 0.5*quadMap_parallel(R_xyu, Hg_mid);

%interval evaluation
for i=1:length(Hf)
    %error_x_rad_old(i,1) = 0.5*dz'*Hf_rad{i}*dz;
    error_x_rad(i,1) = 0.5*dz_abs'*Hf_rad{i}*dz_abs;
end
for i=1:length(Hg)
    %error_y_rad_old(i,1) = 0.5*dz'*Hg_rad{i}*dz;
    error_y_rad(i,1) = 0.5*dz_abs'*Hg_rad{i}*dz_abs;
end



%combine results
error_x_rad_zono = zonotope(interval(-error_x_rad, error_x_rad));
error_y_rad_zono = zonotope(interval(-error_y_rad, error_y_rad));
error_x = error_x_mid + error_x_rad_zono;
error_y = error_y_mid + error_y_rad_zono;

%compute final error: to be CHECKED IF CORRELATION APPLIES
Z_err_x_mid = error_x_mid.Z;
Z_err_x_add_mid = get(obj.linError.CF_inv*error_y_mid,'Z');
error_mid = zonotope(Z_err_x_mid + Z_err_x_add_mid);
error_rad = error_x_rad_zono + obj.linError.CF_inv*error_y_rad_zono;
error = error_mid + error_rad;
%error = error_x + obj.linError.CF_inv*error_y;

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

%------------- END OF CODE --------------