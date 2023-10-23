function derivatives(obj,options)
% derivatives - computes multivariate derivatives (Jacobians, Hessians, etc.)
%    of nonlinear systems in a symbolic way; the result is stored in
%    m-files to which obj has access via handles in its properties
%
% Syntax:
%    derivatives(obj,options)
%
% Inputs:
%    obj - system object
%    options - options struct
%
% Outputs:
%    -
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       29-October-2007 (MA)
% Last update:   07-September-2012 (MA)
%                12-October-2015 (MA)
%                08-April-2016 (MA)
%                10-June-2017 (NK)
%                15-June-2017 (NK)
%                28-June-2017 (NK)
%                06-July-2017
%                15-July-2017 (NK)
%                16-July-2017 (NK)
%                05-November-2017 (MA, generalization for all contDynamics classes)
%                12-November-2017 (MA)
%                03-December-2017 (MA)
%                14-January-2018 (MA)
%                12-November-2018 (NK, removed lagrange remainder files)
%                29-January-2021 (MW, simplify checks, outsource tensor check)
%                18-November-2022 (MW, integrate output equation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set standard path
path = [CORAROOT filesep 'models' filesep 'auxiliary' filesep obj.name];
if ~exist(path,'dir')
   mkdir(path); 
end
addpath(path);

% create symbolic variables
[vars,varsDer] = symVariables(obj,'LRbrackets');

% insert symbolic variables into the system equations
fcon = []; % init constraint equations (non-empty only for nonlinDASys)
if isempty(vars.y) && isempty(vars.p)
    % class: nonlinearSys, nonlinearSysDT
    fdyn = obj.mFile(vars.x,vars.u);
    fout = obj.out_mFile(vars.x,vars.u);
elseif isempty(vars.y) && ~isempty(vars.p)
    % class: nonlinParamSys
    fdyn = obj.mFile(vars.x,vars.u,vars.p);
    fout = obj.out_mFile(vars.x,vars.u,vars.p);
elseif ~isempty(vars.y)
    % class: nonlinDASys
    fdyn = obj.dynFile(vars.x,vars.y,vars.u);
    fcon = obj.conFile(vars.x,vars.y,vars.u);
    fout = obj.out_mFile(vars.x,vars.y,vars.u);
end
    
% check if old derivatives can be re-used
[requiredFiles, requiredFiles_out, allnew, storedData, filepathOld] = ...
    checkTensorRecomputation(obj,fdyn,fcon,fout,vars,options);

% create files:
% 1. allnew = true: remove all files, create new ones
% 2. allnew = false: check which files need to be added

if allnew
    currDir = pwd;
    cd(path);
    delete *.m
    cd(currDir);
end

% already finished if nothing new required
if ~any([allnew;...
        requiredFiles.standard;requiredFiles.int;requiredFiles.higherOrders;...
        requiredFiles_out.standard;requiredFiles_out.int;requiredFiles_out.higherOrders])
    return;
end

fprintf("Computing multivariate derivatives for dynamic '%s':\n", obj.name)

% first, store the actual dynamics (given in symbolic variables)
% and all other Lagrange remainder options in cora/models/auxiliary/<func>
storedData.fdyn = fdyn;
storedData.fcon = fcon;
storedData.fout = fout;

% remove lagrangeRem from storedData
if isfield(storedData,'lagrangeRem')
    storedData = rmfield(storedData,'lagrangeRem');
end

% if lagrangeRem option given, save in storedData
if isfield(options,'lagrangeRem')
    if isfield(options.lagrangeRem,'simplify')
        storedData.lagrangeRem.simplify = options.lagrangeRem.simplify;
    end
    if isfield(options.lagrangeRem,'replacements')
        if ~isempty(vars.p)
            storedData.lagrangeRem.replacements = ...
                options.lagrangeRem.replacements(vars.x,vars.u,vars.p); 
        else
            storedData.lagrangeRem.replacements = ...
                       options.lagrangeRem.replacements(vars.x,vars.u);
        end
    end
    if isfield(options.lagrangeRem,'tensorParallel')
        storedData.lagrangeRem.tensorParallel = ...
                                    options.lagrangeRem.tensorParallel;
    end
    if isfield(options.lagrangeRem,'method')
        storedData.lagrangeRem.method = options.lagrangeRem.method;
    end
end

% nonlinParamSys: save paramInt
if ~isempty(vars.p)
    storedData.paramInt = options.paramInt; 
end
% --- storedData finished

% read value for setting 'simplify' if provided, otherwise choose default
simplify = getDefaultValue('lagrangeRem.simplify',obj,struct(),options,'options');
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'simplify')
    simplify = options.lagrangeRem.simplify;
end

% compute Jacobian
if requiredFiles.standard(1)
    disp('  .. compute symbolic Jacobian');
    name = ['jacobian_' obj.name];
    [Jdyn,Jcon,Jp] = aux_jacobians(obj,fdyn,fcon,vars,options,simplify);
    createJacobianFile(Jdyn,Jcon,Jp,path,name,vars);
    if ~isempty(vars.p)
        if ~isempty(Jp)
            name = ['parametricDynamicFile_' obj.name];
            createParametricDynamicFile(obj,obj.mFile,path,name);
        end
        name = ['jacobian_freeParam_' obj.name];
        createJacobianFile_freeParam(Jdyn,path,name);
    end
end
% compute Jacobian (output)
if requiredFiles_out.standard(1)
    disp('  .. compute symbolic Jacobian (output)');
    name = ['out_jacobian_' obj.name];
    [Jout,~,Jpout] = aux_jacobians(obj,fout,[],vars,options,simplify);
    createJacobianFile(Jout,[],Jpout,path,name,vars);
    if ~isempty(vars.p)
        if ~isempty(Jpout)
            name = ['out_parametricDynamicFile_' obj.name];
            createParametricDynamicFile(obj,obj.out_mFile,path,name);
        end
        name = ['out_jacobian_freeParam_' obj.name];
        createJacobianFile_freeParam(Jout,path,name);
    end
end

% Hessians and third-order tensors
if any([requiredFiles.standard(2:3); requiredFiles.int(2:3)])
    % compute Hessians
    disp('  .. compute symbolic Hessians');
    [J2dyn,J2con] = aux_hessians(fdyn,fcon,vars,simplify);

    if requiredFiles.standard(2)
        % Hessian without interval arithmetic
        name = ['hessianTensor_' obj.name];
        createHessianTensorFile(J2dyn,J2con,path,name,vars,false,options);
    end
    if requiredFiles.int(2)
        % Hessian with interval arithmetic
        name = ['hessianTensorInt_' obj.name];
        createHessianTensorFile(J2dyn,J2con,path,name,vars,true,options);
    end
    
    if any([requiredFiles.standard(3); requiredFiles.int(3)])
        % compute third-order derivatives
        disp('  .. compute symbolic third-order derivatives'); 
        [J3dyn,J3con] = aux_thirdOrderDerivatives(J2dyn,J2con,vars,simplify);
    
        if requiredFiles.standard(3)
            % third-order tensor without interval arithmetic
            name = ['thirdOrderTensor_' obj.name];
            create3rdOrderTensorFile(J3dyn,J3con,path,name,vars,false,options);
        end

        if requiredFiles.int(3)
            % third-order tensor with interval arithmetic
            name = ['thirdOrderTensorInt_' obj.name];
            create3rdOrderTensorFile(J3dyn,J3con,path,name,vars,true,options);
        end
    end

end


% Hessians and third-order tensors (output)
if any([requiredFiles_out.standard(2:3); requiredFiles_out.int(2:3)])
    % compute Hessians
    disp('  .. compute symbolic Hessians (output)');
    J2out = aux_hessians(fout,[],vars,simplify);

    if requiredFiles_out.standard(2)
        % Hessian without interval arithmetic
        name = ['out_hessianTensor_' obj.name];
        createHessianTensorFile(J2out,[],path,name,vars,false,options);
    end
    if requiredFiles_out.int(2)
        % Hessian with interval arithmetic
        name = ['out_hessianTensorInt_' obj.name];
        createHessianTensorFile(J2out,[],path,name,vars,true,options);
    end
    
    if any([requiredFiles_out.standard(3); requiredFiles_out.int(3)])
        % compute third-order derivatives
        disp('  .. compute symbolic third-order derivatives (output)'); 
        J3out = aux_thirdOrderDerivatives(J2out,[],vars,simplify);
    
        if requiredFiles_out.standard(3)
            % third-order tensor without interval arithmetic
            name = ['out_thirdOrderTensor_' obj.name];
            create3rdOrderTensorFile(J3out,[],path,name,vars,false,options);
        end

        if requiredFiles_out.int(3)
            % third-order tensor with interval arithmetic
            name = ['out_thirdOrderTensorInt_' obj.name];
            create3rdOrderTensorFile(J3out,[],path,name,vars,true,options);
        end
    end

end

% 4th and higher-order tensors (not for output equation)
if any(requiredFiles.higherOrders) % equivalent to options.tensorOrder >= 4
    createHigherOrderTensorFiles(fdyn,vars,varsDer,path,obj.name,options); 
end

% rehash the folder so that new generated files are used
rehash path;

% save data so that symbolic computations do not have to be re-computed
save(filepathOld,'storedData');

% 
fprintf("Done.\n")

end


% Auxiliary functions -----------------------------------------------------

function [Jdyn,Jcon,Jp] = aux_jacobians(obj,fdyn,fcon,vars,options,simplifyOpt)
% jacobians - compute symbolic Jacobians of differential equation,
%    constraint equation and w.r.t parameter uncertainties
%
% Inputs:
%    obj - contDynamics object
%    fdyn - symbolic differential equation
%    fcon - symbolic constraint equation (only nonlinDASys)
%    vars - symbolic variables
%    options - options for reachability analysis, mainly
%                options.tensorOrder and options.lagrangeRem
%    simplifyOpt - specification if and how tensor should be simplified
%
% Outputs:
%    Jdyn - Jacobian of differential equation
%    Jcon - Jacobian of constraint equation (only nonlinDASys)
%    Jp - Jacobian of parameters (only nonlinParamSys)

% init
Jdyn = [];
Jcon = [];
Jp = [];

%compute jacobian with respect to the state
% dynamic part
if ~isempty(fdyn)
    Jdyn.x = jacobian(fdyn,vars.x);
else
    Jdyn.x = [];
end
% constraint part
if ~isempty(fcon)
    Jcon.x = jacobian(fcon,vars.x);
end

%compute jacobian with respect to the input
% dynamic part
if ~isempty(fdyn)
    Jdyn.u = jacobian(fdyn,vars.u);
else
    Jdyn.u = [];
end
% constraint part
if ~isempty(fcon)
    Jcon.u = jacobian(fcon,vars.u);
end

%compute jacobian with respect to the constraint state
if ~isempty(fcon)
    % dynamic part
    Jdyn.y = jacobian(fdyn,vars.y);
    % constraint part
    Jcon.y = jacobian(fcon,vars.y);
end

% perform simplification
if strcmp(simplifyOpt,'simplify') 
    Jdyn.x = simplify(Jdyn.x);
    Jdyn.u = simplify(Jdyn.u);
    if ~isempty(fcon)
        Jcon.x = simplify(Jcon.x);
        Jdyn.y = simplify(Jdyn.y);
        Jcon.y = simplify(Jcon.y);
        Jcon.u = simplify(Jcon.u);
    end
elseif strcmp(simplifyOpt,'collect')
    Jdyn.x = collect(Jdyn.x,vars.x);
    Jdyn.u = collect(Jdyn.u,vars.x);
    if ~isempty(fcon)
        Jcon.x = collect(Jcon.x);
        Jdyn.y = collect(Jdyn.y);
        Jcon.y = collect(Jcon.y);
        Jcon.u = collect(Jcon.u);
    end
end

% special derivatives for nonlinear systems with parameters linearly
% influencing the derivative
if ~isempty(vars.p)
    
    %store jacobians with respect to parameters
    Jp.x=cell(1,obj.nrOfParam+1);
    Jp.u=cell(1,obj.nrOfParam+1);

    %part without parameters
    try
        Jp.x{1} = subs(Jdyn.x,vars.p,zeros(obj.nrOfParam,1));
        Jp.u{1} = subs(Jdyn.u,vars.p,zeros(obj.nrOfParam,1));
        %part with parameters
        I = eye(obj.nrOfParam); %identity matrix
        for i=1:obj.nrOfParam
            Jp.x{i+1} = subs(Jdyn.x,vars.p,I(:,i)) - Jp.x{1};
            Jp.u{i+1} = subs(Jdyn.u,vars.p,I(:,i)) - Jp.u{1};
        end

        % if parameters are uncertain within an interval
        if isa(options.paramInt,'interval')
            %normalize
            pCenter = center(options.paramInt);
            pDelta = rad(options.paramInt);

            for i=1:obj.nrOfParam
                %center
                Jp.x{1} = Jp.x{1} + pCenter(i)*Jp.x{i+1};
                Jp.u{1} = Jp.u{1} + pCenter(i)*Jp.u{i+1};
                %generators
                Jp.x{i+1} = pDelta(i)*Jp.x{i+1};
                Jp.u{i+1} = pDelta(i)*Jp.u{i+1};
            end
        else
            for i=1:obj.nrOfParam
                %center
                Jp.x{1} = Jp.x{1} + options.paramInt(i)*Jp.x{i+1};
                Jp.u{1} = Jp.u{1} + options.paramInt(i)*Jp.u{i+1};
            end
        end
    catch
        Jp = [];
        disp('Parameters are not linearly influencing the system.');
    end
end

end

function [Hdyn, Hcon] = aux_hessians(fdyn,fcon,vars,simplifyOpt)
% hessians - compute Hessian tensors for differential equation and
%    constraint equation (only nonlinDASys)
%
% Inputs:
%    fdyn - symbolic differential equation
%    fcon - symbolic constraint equation (only nonlinDASys)
%    vars - symbolic variables
%    simplifyOpt - specification if and how tensor should be simplified
%
% Outputs:
%    Hdyn - Hessian tensor of differential equation
%    Hcon - Hessian tensor of constraint equation (only nonlinDASys)

% init
Hdyn = sym([]);
Hcon = sym([]);

%compute second-order jacobians using 'LR' variables
if isempty(vars.y) % no constraint equations
    vars.z = [vars.x;vars.u];
else % with constraint equations
    vars.z = [vars.x;vars.y;vars.u];
end

%dynamic jacobians
Jdyn_comb = jacobian(fdyn,vars.z);
for k=1:length(Jdyn_comb(:,1))
    % compute 2nd- order Jacobians
    Hdyn(k,:,:)=jacobian(Jdyn_comb(k,:),vars.z);
end

%constraint jacobians
if ~isempty(fcon)
    Jcon_comb = jacobian(fcon,vars.z);
    for k=1:length(Jcon_comb(:,1))
        %Calculate 2nd order Jacobians
        Hcon(k,:,:)=jacobian(Jcon_comb(k,:),vars.z);
    end
end

%potentially simplify expression         
if strcmp(simplifyOpt,'simplify') 
    for k=1:length(Jdyn_comb(:,1))
        Hdyn(k,:,:) = simplify(Hdyn(k,:,:));
        if ~isempty(Hcon)
            Hcon(k,:,:) = simplify(Hcon(k,:,:));
        end
    end
elseif strcmp(simplifyOpt,'collect')
    for k=1:length(Jdyn_comb(:,1))
        Hdyn(k,:,:) = collect(Hdyn(k,:,:),vars.x);
        if ~isempty(Hcon)
            Hcon(k,:,:) = collect(Hcon(k,:,:),vars.x);
        end
    end
end
    
end

function [J3dyn, J3con] = aux_thirdOrderDerivatives(J2dyn,J2con,vars,simplifyOpt)
% thirdOrderDerivatives - compute third-order derivatives for
%    differential equation and constraint equation (only nonlinDASys)
%
% Inputs:
%    J2dyn - symbolic Hessian tensor of differential equation
%    J2con - symbolic Hessian tensor of constraint equation (only nonlinDASys)
%    vars - symbolic variables
%    simplifyOpt - specification if and how tensor should be simplified
%
% Outputs:
%    J3dyn - symbolic third-order tensor of differential equation
%    J3con - symbolic third-order tensor of constraint equation (only nonlinDASys)


dim = length(J2dyn(:,1,1));
nrOfVars = length(J2dyn(1,:,1));
J3dyn = sym(zeros(dim,nrOfVars,nrOfVars,nrOfVars));
J3con = sym(zeros(dim,nrOfVars,nrOfVars,nrOfVars));

% construct vector for which derivative is computed
if isempty(vars.y) % no constraint equations
    vars.z = [vars.x;vars.u];
else % with constraint equations
    vars.z = [vars.x;vars.y;vars.u];
end

%compute third-order jacobians using 'LR' variables
% dynamic part
for k=1:length(J2dyn(:,1,1))
    for l=1:length(J2dyn(1,:,1))
        % compute 3rd-order Jacobians
        if ~isempty(find(J2dyn(k,l,:), 1))
            J3dyn(k,l,:,:)=jacobian(reshape(J2dyn(k,l,:),[nrOfVars,1]),vars.z);
        end
    end
end
% constraint part
if ~isempty(J2con)
    for k=1:length(J2con(:,1,1))
        for l=1:length(J2con(1,:,1))
            % compute 3rd-order Jacobians
            if ~isempty(find(J2con(k,l,:), 1))
                J3con(k,l,:,:)=jacobian(reshape(J2con(k,l,:),[nrOfVars,1]),vars.z);
            end
        end
    end
end

% potentially simplify expression            
if strcmp(simplifyOpt,'simplify') 
    for k=1:length(J2dyn(:,1,1))
        for l=1:length(J2dyn(1,:,1))
            if ~isempty(find(J2dyn(k,l,:), 1))
                J3dyn(k,l,:,:) = simplify(J3dyn(k,l,:,:));
                if ~isempty(J3con)
                    J3con(k,l,:,:) = simplify(J3con(k,l,:,:));
                end
            end
        end
    end
elseif strcmp(simplifyOpt,'collect')
    for k=1:length(J2dyn(:,1,1))
        for l=1:length(J2dyn(1,:,1))
            if ~isempty(find(J2dyn(k,l,:), 1))
                J3dyn(k,l,:,:) = collect(J3dyn(k,l,:,:),vars.x);
                if ~isempty(J3con)
                    J3con(k,l,:,:) = collect(J3con(k,l,:,:),vars.x);
                end
            end
        end
    end
end
    
end

% ------------------------------ END OF CODE ------------------------------
