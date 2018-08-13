function [errorTotal, errorInt, errorInt_x, errorInt_y, Rtotal_y] = linError_thirdOrder_comp(obj, options, R, Verror_y)
% linError_thirdOrder - computes the linearization error using a third
% order Taylor expansion
%
% Syntax:  
%    [error, errorInt, errorInt_x, errorInt_y, R_y] = linError_thirdOrder(obj, options, R, Verror_y)
%
% Inputs:
%    obj - nonlinear DAE system object
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
% Written:      26-June-2013
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%compute set of algebraic variables
f0_con = obj.linError.f0_con;
D = obj.linError.D;
E = obj.linError.E;
F_inv = obj.linError.F_inv;
R_y_cor = -F_inv*(f0_con + D*R); %correlated part
R_y_add = -F_inv*(E*options.U + Verror_y); %uncorrelated part
R_y = R_y_cor + R_y_add;

%compute variable projections onto subsystems
if strcmp(options.mode,'normal')
    subsystem = options.subsystem;
elseif strcmp(options.mode,'fault')
    subsystem = options.subsystemFault;
end

%init
errorTotal = zeros(obj.dim, 1);
errorTotal_x = zeros(obj.dim, 1);
errorTotal_y = zeros(obj.nrOfConstraints, 1);

%obtain indices for projection onto subsystems
index = indexForSubsystems(obj, subsystem, 3);

%linearization errors
for iSys = 1:length(index)
    psub{iSys}.x = index{iSys}.X * obj.linError.p.x;
    psub{iSys}.y = index{iSys}.Y * obj.linError.p.y;
    psub{iSys}.u = index{iSys}.Uv*obj.other.Vgen' + index{iSys}.Uy*obj.linError.p.y + index{iSys}.Uu*obj.linError.p.u; %IMPORTANT: INCLUDE effect of Vgen
end

for iSys = 1:length(index)
    %project
    %reachable sets
    Rsub = index{iSys}.X * R;
    Rsub_y_cor = index{iSys}.Y * R_y_cor;
    Rsub_y_add = index{iSys}.Y * R_y_add;
    Usub = index{iSys}.Uy*R_y + index{iSys}.Uu*options.U; %IMPORTANT: do not include effect of Vgen
    
    %obtain intervals and combined interval z
    dx = interval(Rsub);
    dy = interval(Rsub_y_cor + Rsub_y_add);
    du = interval(Usub);
    dz = [dx; dy; du];
    
    %compute interval of reachable set
    totalInt_x = dx + psub{iSys}.x;

    %compute interval of algebraic states
    totalInt_y = dy + psub{iSys}.y;

    %compute intervals of input
    totalInt_u = du + psub{iSys}.u;

    %obtain hessian and third order tensor
    [Hf, Hg] = subsystem{iSys}.hessian(obj.linError.p.x, obj.linError.p.y, obj.linError.p.u);
    [Tf, Tg] = subsystem{iSys}.thirdOrderTensor(totalInt_x, totalInt_y, totalInt_u);
    
    %compute zonotope of state, constarint variables, and input
    Z_x = get(Rsub,'Z');
    Z_y_cor = get(Rsub_y_cor,'Z');
    Z_y_add = get(Rsub_y_add,'Z');
    Z_0 = zeros(length(Z_x(:,1)), length(Z_y_add(1,:)));
    R_xy = zonotope([Z_x, Z_0; Z_y_cor, Z_y_add]);
    R_xyu = cartesianProduct(R_xy, Usub);
    R_xyu = reduce(R_xyu,'girard',options.errorOrder);

    %obtain absolute values
    dz_abs = max(abs(infimum(dz)), abs(supremum(dz)));
    
    %separate evaluation
    Hf_real = []; %delete previous values
    Hg_real = []; %delete previous values
    
    %store Hf and Hg as real-valued 
    for i=1:length(Hf)
        Hf_real{i} = mid(Hf{i});
    end
    for i=1:length(Hg)
        Hg_real{i} = mid(Hg{i});
    end
    
    %zonotope evaluation
    %algebraic part
    error_y_secondOrder = 0.5*quadraticMultiplication_parallel(R_xyu, Hg_real);
    %third order interval evaluation (algebraic part)
    for i=1:length(Tg(:,1))
        error_sum = interval(0,0);
        for j=1:length(Tg(1,:))
            error_tmp = dz'*Tg{i,j}*dz;
            error_sum = error_sum + error_tmp * dz(j);
        end
        error_y_thirdOrder(i,1) = 1/6*error_sum;
    end

    %combine results
    error_thirdOrder_y_zono = zonotope(error_y_thirdOrder);
    error_y{iSys} = error_y_secondOrder + error_thirdOrder_y_zono;
    
    %there have to be state variables; otherwise the error for x is 0
    if ~isempty(Hf_real)
        error_x_secondOrder = 0.5*quadraticMultiplication_parallel(R_xyu, Hf_real);
        error_x_thirdOrder = intval();
        %third order interval evaluation (dynamic part)
        for i=1:length(Tf(:,1))
            error_sum = interval(0,0);
            for j=1:length(Tf(1,:))
                error_tmp = dz'*Tf{i,j}*dz;
                error_sum = error_sum + error_tmp * dz(j);
            end
            error_x_thirdOrder(i,1) = 1/6*error_sum;
        end
        
        %combine results
        error_thirdOrder_x_zono = zonotope(error_x_thirdOrder);
        error_x{iSys} = error_x_secondOrder + error_thirdOrder_x_zono;
        
        %compute auxiliary values
        [Asub,Bsub,Csub,Dsub,Esub,Fsub] = subsystem{iSys}.jacobian(psub{iSys}.x, psub{iSys}.y, psub{iSys}.u);
        Fsub_inv = pinv(Fsub);
        CFsub_inv = Csub*Fsub_inv;
        
        %compute final error
        Z_err_x = get(error_x_secondOrder,'Z');
        Z_err_x_add = get(CFsub_inv*error_y_secondOrder,'Z');
        error_secondOrder = zonotope(Z_err_x + Z_err_x_add);
        error_thirdOrder = error_thirdOrder_x_zono + CFsub_inv*error_thirdOrder_y_zono;
        error{iSys} = error_secondOrder + error_thirdOrder;
    end
end

%compute total errors
for iSys = 1:length(index)
    if ~isempty(index{iSys}.X)
        %reduce
        error{iSys} = reduce(error{iSys},'girard',options.zonotopeOrder);
        error_x{iSys} = reduce(error_x{iSys},'girard',options.zonotopeOrder);
        error_y{iSys} = reduce(error_y{iSys},'girard',options.zonotopeOrder);
        
        errorTotal = errorTotal + index{iSys}.X'*error{iSys};
        errorTotal_x = errorTotal_x + index{iSys}.X'*error_x{iSys};
    end
    errorTotal_y = errorTotal_y + index{iSys}.Y'*error_y{iSys};
end

%reduce
errorTotal = reduce(errorTotal,'girard',options.zonotopeOrder);
errorTotal_y = reduce(errorTotal_y,'girard',options.zonotopeOrder);

%update R_y
Rtotal_y =  obj.linError.p.y + (-F_inv)*(f0_con + D*R + E*options.U + errorTotal_y);

%error intervals
errorIHabs = abs(interval(errorTotal));
errorInt = supremum(errorIHabs);

errorIHabs_y = abs(interval(errorTotal_y));
errorInt_y = supremum(errorIHabs_y);

errorIHabs_x = abs(interval(errorTotal_x));
errorInt_x = supremum(errorIHabs_x);


%------------- END OF CODE --------------