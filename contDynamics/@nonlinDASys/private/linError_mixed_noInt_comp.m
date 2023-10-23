function [errorTotal, errorInt, errorInt_x, errorInt_y, Rtotal_y] = ...
    linError_mixed_noInt_comp(obj, options, R, Verror_y)
% linError_mixed_noInt_comp - computes the linearization error 
% compositionally according to [1].
%
% Syntax:
%    [errorTotal, errorInt, errorInt_x, errorInt_y, Rtotal_y] = ...
%               linError_mixed_noInt_comp(obj, options, R, Verror_y)
%
% Inputs:
%    obj - nonlinear differential algebraic system object
%    options - options struct
%    R - actual reachable set
%    Verror_y - set of algebraic linearization error
%
% Outputs:
%    errorTotal - zonotope overapproximating the linearization error
%    errorInt - interval overapproximating the linearization error
%    errorInt_x - interval overapproximating the linearization error (dynamic part)
%    errorInt_y - interval overapproximating the linearization error (constraint part)
%    Rtotal_y - reachable set of the algebraic part
%
% Reference:
%    [1] M. Althoff, "Formal and Compositional Analysis of Power Systems 
%        using Reachable Sets", IEEE Transactions on Power Systems 29 (5), 
%        2014, 2270-2280
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       21-November-2011
% Last update:   23-May-2013
%                06-June-2013
%                04-October-2015
%                25-July-2016 (intervalhull replaced by interval)
%                08-June-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute set of algebraic variables
f0_con = obj.linError.f0_con;
D = obj.linError.D;
E = obj.linError.E;
F_inv = obj.linError.F_inv;
R_y_cor = -F_inv*(f0_con + D*R); %correlated part
R_y_add = -F_inv*(E*options.U + Verror_y); %uncorrelated part
R_y = R_y_cor + R_y_add;

% init
errorTotal = zeros(obj.dim, 1);
errorTotal_x = zeros(obj.dim, 1);
errorTotal_y = zeros(obj.nrOfConstraints, 1);

% obtain indices for projection onto subsystems
index = options.index;

% linearization errors
for iSys = 1:length(index)
    psub{iSys}.x = index{iSys}.X * obj.linError.p.x;
    psub{iSys}.y = index{iSys}.Y * obj.linError.p.y;
    psub{iSys}.u = index{iSys}.Uv*options.Vgen + index{iSys}.Uy*obj.linError.p.y + index{iSys}.Uu*obj.linError.p.u; %IMPORTANT: INCLUDE effect of Vgen
end

parfor iSys = 1:length(index)
%for iSys = 1:length(index)
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

    %obtain hessian tensor
    [Hf, Hg] = options.subsystemHessian{iSys}(totalInt_x, totalInt_y, totalInt_u);

    %compute zonotope of state, constarint variables, and input
    Z_x = [Rsub.c,Rsub.G];
    Z_y_cor = [Rsub_y_cor.c,Rsub_y_cor.G];
    Z_y_add = [Rsub_y_add.c,Rsub_y_add.G];
    Z_0 = zeros(length(Z_x(:,1)), length(Z_y_add(1,:)));
    R_xy = zonotope([Z_x, Z_0; Z_y_cor, Z_y_add]);
    R_xyu = cartProd(R_xy, Usub);
    R_xyu = reduce(R_xyu,'girard',options.errorOrder);


    %obtain absolute values
    %dz_abs = max(abs(inf(dz)), abs(sup(dz)));
    dz_abs = max(abs(infimum(dz)), abs(supremum(dz)));

    %separate evaluation
    Hf_mid = []; %delete previous values
    Hg_mid = []; %delete previous values
    Hf_rad = []; %delete previous values
    Hg_rad = []; %delete previous values
    for i=1:length(Hf)
        Hf_mid{i} = sparse(center(Hf{i}));
        Hf_rad{i} = sparse(rad(Hf{i}));
    end
    for i=1:length(Hg)
        Hg_mid{i} = sparse(center(Hg{i}));
        Hg_rad{i} = sparse(rad(Hg{i}));
    end
    %zonotope evaluation
    %algebraic part
    error_y_mid = 0.5*quadMap_parallel(R_xyu, Hg_mid);
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
        error_x_mid = 0.5*quadMap_parallel(R_xyu, Hf_mid);
        %interval evaluation
        error_x_rad = []; %delete previous values
        for i=1:length(Hf)
            error_x_rad(i,1) = 0.5*dz_abs'*Hf_rad{i}*dz_abs;
        end
        %combine results
        error_x_rad_zono = zonotope(interval(-error_x_rad, error_x_rad));
        error_x{iSys} = error_x_mid + error_x_rad_zono;
        
        %compute auxiliary values
        [~,~,Csub,~,~,Fsub] = options.subsystemJacobian{iSys}(psub{iSys}.x, psub{iSys}.y, psub{iSys}.u);
        Fsub_inv = pinv(Fsub);
        CFsub_inv = Csub*Fsub_inv;
        
        %compute final error
        Z_err_x_mid = [error_x_mid.c,error_x_mid.G];
        Z_err_x_add_mid = [CFsub_inv*error_y_mid.c,CFsub_inv*error_y_mid.G];
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

% ------------------------------ END OF CODE ------------------------------
