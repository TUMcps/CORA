function [errorTotal, errorInt, errorInt_x, errorInt_y, Rtotal_y] = linError_mixed_noInt_comp(obj, options, R, Verror_y)
% linError_mixed_noInt_comp - computes the linearization error compositionally
%
% Syntax:  
%    [obj] = linError(obj,options)
%
% Inputs:
%    obj - nonlinear DAE system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    error - zonotope overapproximating the linearization error
%    errorInt - interval overapproximating the linearization error
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
%               06-June-2013
%               15-June-2016
%               25-July-2016 (intervalhull replaced by interval)
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

%parfor iSys = 1:length(index)
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

    %obtain absolute values
    dx_abs = max(abs(infimum(dx)), abs(supremum(dx)));
    dy_abs = max(abs(infimum(dy)), abs(supremum(dy)));
    du_abs = max(abs(infimum(du)), abs(supremum(du)));
    dz_abs = max(abs(infimum(dz)), abs(supremum(dz)));

    %obtain hessian tensor
    [Hf, Hg] = subsystem{iSys}.hessian(psub{iSys}.x, psub{iSys}.y, psub{iSys}.u);

    %compute zonotope of state, constarint variables, and input
    Z_x = get(Rsub,'Z');
    Z_y_cor = get(Rsub_y_cor,'Z');
    Z_y_add = get(Rsub_y_add,'Z');
    Z_0 = zeros(length(Z_x(:,1)), length(Z_y_add(1,:)));
    R_xy = zonotope([Z_x, Z_0; Z_y_cor, Z_y_add]);
    R_xyu = cartesianProduct(R_xy, Usub);
    R_xyu = reduce(R_xyu,'girard',options.errorOrder);



    %separate evaluation
    Hf_mid = []; %delete previous values
    Hg_mid = []; %delete previous values
    Hf_rad = []; %delete previous values
    Hg_rad = []; %delete previous values
    for i=1:length(Hf)
        Hf_mid{i} = sparse(mid(Hf{i}));
        Hf_rad{i} = sparse(rad(Hf{i}));
    end
    for i=1:length(Hg)
        Hg_mid{i} = sparse(mid(Hg{i}));
        Hg_rad{i} = sparse(rad(Hg{i}));
    end
    %zonotope evaluation
    %algebraic part
    error_y_mid = 0.5*quadraticMultiplication_parallel(R_xyu, Hg_mid);
    %interval evaluation
    error_y_rad = []; %delete previous values
    for i=1:length(Hg)
        error_y_rad(i,1) = 0.5*dz_abs'*Hg_rad{i}*dz_abs;
    end

    %combine results
    error_y_rad_zono = zonotope(interval(-error_y_rad, error_y_rad));
    error_y{iSys} = error_y_mid + error_y_rad_zono;
    
    %there have to be state variables; otherwise the error for x is 0
    if ~isempty(Hf_mid)
        error_x_mid = 0.5*quadraticMultiplication_parallel(R_xyu, Hf_mid);
        %interval evaluation
        error_x_rad = []; %delete previous values
        for i=1:length(Hf)
            error_x_rad(i,1) = 0.5*dz_abs'*Hf_rad{i}*dz_abs;
        end
        %combine results
        error_x_rad_zono = zonotope(interval(-error_x_rad, error_x_rad));
        error_x{iSys} = error_x_mid + error_x_rad_zono;
        
        %compute auxiliary values
        [Asub,Bsub,Csub,Dsub,Esub,Fsub] = subsystem{iSys}.jacobian(psub{iSys}.x, psub{iSys}.y, psub{iSys}.u);
        Fsub_inv = pinv(Fsub);
        CFsub_inv = Csub*Fsub_inv;
        
        %compute final error
        Z_err_x_mid = get(error_x_mid,'Z');
        Z_err_x_add_mid = get(CFsub_inv*error_y_mid,'Z');
        error_mid = zonotope(Z_err_x_mid + Z_err_x_add_mid);
        error_rad = error_x_rad_zono + CFsub_inv*error_y_rad_zono;
        error{iSys} = error_mid + error_rad;
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