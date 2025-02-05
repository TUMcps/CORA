function [error,errorInt,errorInt_x,errorInt_y,R_y] = priv_linError_mixed_noInt(nlnsysDA,R,Verror_y,params,options)
% priv_linError_mixed_noInt - computes the linearization error
%
% Syntax:
%    [error,errorInt,errorInt_x,errorInt_y,R_y] = ...
%           priv_linError_mixed_noInt(nlnsysDA,options,R,Verror_y)
%
% Inputs:
%    nlnsysDA - nonlinDASys object
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
% Written:       21-November-2011
% Last update:   23-May-2013
%                25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

%obtain hessian tensor
nlnsysDA.setHessian('int');
[Hf, Hg] = nlnsysDA.hessian(totalInt_x, totalInt_y, totalInt_u);

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
    error_x_rad(i,1) = 0.5*dz_abs'*Hf_rad{i}*dz_abs;
end
for i=1:length(Hg)
    error_y_rad(i,1) = 0.5*dz_abs'*Hg_rad{i}*dz_abs;
end

%combine results
error_x_rad_zono = zonotope(interval(-error_x_rad, error_x_rad));
error_y_rad_zono = zonotope(interval(-error_y_rad, error_y_rad));
error_x = error_x_mid + error_x_rad_zono;
error_y = error_y_mid + error_y_rad_zono;

%compute final error: to be CHECKED IF CORRELATION APPLIES
Z_err_x_mid = [error_x_mid.c,error_x_mid.G];
Z_err_x_add_mid = nlnsysDA.linError.CF_inv*[error_y_mid.c,error_y_mid.G];
error_mid = zonotope(Z_err_x_mid + Z_err_x_add_mid);
error_rad = error_x_rad_zono + nlnsysDA.linError.CF_inv*error_y_rad_zono;
error = error_mid + error_rad;
%error = error_x + obj.linError.CF_inv*error_y;

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
