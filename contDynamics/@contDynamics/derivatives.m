function derivatives(varargin)
% derivatives - computes multivariate derivatives (jacobians, hessians, etc.)
% of nonlinear systems in a symbolic way; the result is stored in m-files 
% and passed by a handle
%
% Syntax:  
%    derivatives(varargin)
%
% Inputs:
%    obj - system object
%    options - options struct
%
% Outputs:
%    -
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      29-October-2007 (MA)
% Last update:  07-September-2012 (MA)
%               12-October-2015 (MA)
%               08-April-2016 (MA)
%               10-June-2017 (NK)
%               15-June-2017 (NK)
%               28-June-2017 (NK)
%               06-July-2017
%               15-July-2017 (NK)
%               16-July-2017 (NK)
%               05-November-2017 (MA, generalization for all contDynamics classes)
%               12-November-2017 (MA)
%               03-December-2017 (MA)
%               14-January-2018 (MA)
%               12-November-2018 (NK, removed lagrange remainder files)
% Last revision:---

%------------- BEGIN CODE --------------

% obtain object and options
obj = varargin{1};
options = varargin{2};

% set standard path
path = [coraroot '/models/auxiliary/',obj.name];
if ~exist(path,'dir')
   mkdir(path); 
end
addpath(path);

% create symbolic variables
[vars,varsDer] = symVariables(obj,'LRbrackets');

% insert symbolic variables into the system equations
fcon = []; %constraint equations (only for differential algebraic equations)
if isempty(vars.y) && isempty(vars.p)
    fdyn = obj.mFile(vars.x,vars.u);
elseif isempty(vars.y) && ~isempty(vars.p)
    fdyn = obj.mFile(vars.x,vars.u,vars.p);
elseif ~isempty(vars.y)
    fdyn = obj.dynFile(vars.x,vars.y,vars.u);
    fcon = obj.conFile(vars.x,vars.y,vars.u);
end
    
% check if old derivations can be used
[generateFiles, storedData, filepathOld] = checkForNewComputation(obj,fdyn,fcon,vars,options);

% create files
if generateFiles
    
    % store the actual dynamics (given in symbolic variables) and the 
    % values for all other lagrange remainder options in the 
    % /StateSpaceModel-directory inside cora
    storedData.fdyn = fdyn;
    storedData.fcon = fcon;

    storedData.tensorOrder = options.tensorOrder;

    if isfield(storedData,'lagrangeRem')
        storedData = rmfield(storedData,'lagrangeRem');
    end
    
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
    
    % get value for setting 'simplify'
    simplify = 'none';
    if isfield(options,'lagrangeRem')
       if isfield(options.lagrangeRem,'simplify')
          simplify = options.lagrangeRem.simplify; 
       end
    end   
    
    % compute jacobians
    disp('create Jacobians');
    [Jdyn,Jcon,Jp] = jacobians(fdyn,fcon,vars,obj,options,simplify);

    % compute hessians
    disp('create hessians');
    [J2dyn, J2con] = hessians(fdyn,fcon,vars,obj,simplify);

    % compute third order derivatives
    if options.tensorOrder>2
        disp('create 3rd order derivatives'); 
        [J3dyn, J3con] = thirdOrderDerivatives(J2dyn,J2con,vars,simplify);     
    end
    
    % generate mFile that computes the Lagrange remainder and the jacobian
    disp('create mFiles')

    createJacobianFile(Jdyn,Jcon,Jp,path,obj.name,vars);
    if ~isempty(vars.p)
        if ~isempty(Jp)
            createParametricDynamicFile(obj,path,obj.name);
        end
        createJacobianFile_freeParam(Jdyn,path,obj.name);
    end
    
    if options.tensorOrder>2
        createHessianTensorFile(J2dyn,J2con,path,obj.name,vars,0,options);
        create3rdOrderTensorFile(J3dyn,J3con,path,obj.name,vars,options);
    else
        createHessianTensorFile(J2dyn,J2con,path,obj.name,vars,1,options);
    end
    
    if options.tensorOrder>3
        createHigherOrderTensorFiles(fdyn,vars,varsDer,path,obj.name,options); 
    end

    % rehash the folder so that new generated files are used
    rehash path;

    % save data so that symbolic computations do not have to be re-computed
    save(filepathOld,'storedData');
end
end

% functions that checks whether symbolic computatios have to be performed
% or weather the old derivations have remained unchanged
function [generateFiles, storedData, filepathOld] = checkForNewComputation(obj,fdyn,fcon,vars,options)

% set standard path
path = [coraroot '/models/auxiliary/' obj.name];

% check if dynamics files already exists
filepathOld = [path filesep obj.name '_lastVersion.mat'];

% init
generateFiles = 1;
storedData = [];

if exist(filepathOld,'file')
    % load stored dynamics 
    load(filepathOld);
    
    % compare stored and actual dynamics (given in symbolic variables)
    if exist('storedData','var')
        equalDynEquations = isequal(fdyn,storedData.fdyn);
        equalConEquations = isequal(fcon,storedData.fcon);
        equalTensorOrder = (options.tensorOrder == storedData.tensorOrder);
        
        % check if equations and tensor order match
        if equalDynEquations && equalConEquations && equalTensorOrder
            
            generateFiles = 0;
            
            if isfield(storedData,'lagrangeRem')
                
                if isfield(options,'lagrangeRem')
                   
                    temp1 = options.lagrangeRem;
                    temp2 = storedData.lagrangeRem;
                    
                    % check setting 'simplify'
                    if isfield(temp1,'simplify')
                        if isfield(temp2,'simplify')
                            if ~strcmp(temp1.simplify,temp2.simplify)
                                generateFiles = 1;
                                return;
                            end
                        else
                            if ~strcmp(temp1.simplify,'none')
                                generateFiles = 1;
                                return; 
                            end
                        end
                    else
                        if isfield(temp2,'simplify') && ~strcmp(temp2.simplify,'none')
                            generateFiles = 1;
                            return;
                        end
                    end
                    
                    % check setting 'tensorParallel'
                    if isfield(temp1,'tensorParallel')
                        if isfield(temp2,'tensorParallel')
                            if temp1.tensorParallel ~= temp2.tensorParallel
                                generateFiles = 1;
                                return; 
                            end
                        else
                            if temp1.tensorParallel ~= 0
                                generateFiles = 1;
                                return; 
                            end
                        end
                    else
                        if isfield(temp2,'tensorParallel') && temp2.tensorParallel ~= 0
                            generateFiles = 1;
                            return;
                        end
                    end
                    
                    % check setting 'method'
                    if isfield(temp1,'method')
                        if isfield(temp2,'method')
                            if ~strcmp(temp1.method,temp2.method)
                                generateFiles = 1;
                                return;
                            end
                        else
                            if ~strcmp(temp1.method,'interval')
                                generateFiles = 1;
                                return; 
                            end
                        end
                    else
                        if isfield(temp2,'method') && ~strcmp(temp2.method,'interval')
                            generateFiles = 1;
                            return;
                        end
                    end
                    
                    % check setting 'replacements'
                    if isfield(temp1,'replacements')
                        if isfield(temp2,'replacements')
                            if ~isempty(vars.p)
                                if ~isequal(temp1.replacements(vars.x,vars.u,vars.p),temp2.replacements)
                                    generateFiles = 1;
                                    return;
                                end
                            else
                                if ~isequal(temp1.replacements(vars.x,vars.u),temp2.replacements)
                                    generateFiles = 1;
                                    return;
                                end
                            end
                        else
                            generateFiles = 1;
                            return;
                        end
                    else
                        if isfield(temp2,'replacements')
                            generateFiles = 1;
                            return;
                        end
                    end
                   
                else
                    
                    % check if options differ from the default values
                    temp = storedData.lagrangeRem;
                    
                    if (isfield(temp,'simplify') && ~strcmp(temp.simplify,'none')) || ...
                       (isfield(temp,'tensorParallel') && ~(temp.tensorParallel == 1)) || ...
                       (isfield(temp,'method') && ~strcmp(temp.method,'interval')) || ...
                       isfield(temp,'replacements')
                        generateFiles = 1;
                    end
                end 
                
            else
                if isfield(options,'lagrangeRem')
                    
                    % check if options differ from the default values
                    temp = options.lagrangeRem;
                    
                    if (isfield(temp,'simplify') && ~strcmp(temp.simplify,'none')) || ...
                       (isfield(temp,'tensorParallel') && ~(temp.tensorParallel == 1)) || ...
                       (isfield(temp,'method') && ~strcmp(temp.method,'interval')) || ...
                       isfield(temp,'replacements')
                        generateFiles = 1;
                    end
                end 
            end
        end
    end
end
end


% compute jacobians
function [Jdyn,Jcon,Jp] = jacobians(fdyn,fcon,vars,obj,options,simplify)

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
if strcmp(simplify,'simplify') 
    Jdyn.x = simplify(Jdyn.x);
    Jdyn.u = simplify(Jdyn.u);
    if ~isempty(fcon)
        Jcon.x = simplify(Jcon.x);
        Jdyn.y = simplify(Jdyn.y);
        Jcon.y = simplify(Jcon.y);
        Jcon.u = simplify(Jcon.u);
    end
elseif strcmp(simplify,'collect')
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

    %store jacobians with respect to parameters (assumption: expressions are linear in parameters)
    %init
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
        disp('parameters are not linearly influencing the system');
    end
end
end

% compute hessians
function [J2dyn, J2con] = hessians(fdyn,fcon,vars,obj,simplify)

    % init
    J2dyn = sym([]);
    J2con = sym([]);

    %compute second order jacobians using 'LR' variables
    if isempty(fcon) % no constraint equations
        vars.z = [vars.x;vars.u];
    else % with constraint equations
        vars.z = [vars.x;vars.y;vars.u];
    end

    %dynamic jacobians
    Jdyn_comb = jacobian(fdyn,vars.z);
    for k=1:length(Jdyn_comb(:,1))
        %Calculate 2nd order Jacobians
        J2dyn(k,:,:)=jacobian(Jdyn_comb(k,:),vars.z);
    end

    %constraint jacobians
    if ~isempty(fcon)
        Jcon_comb = jacobian(fcon,vars.z);
        for k=1:length(Jcon_comb(:,1))
            %Calculate 2nd order Jacobians
            J2con(k,:,:)=jacobian(Jcon_comb(k,:),vars.z);
        end
    end

    %potentially simplify expression         
    if strcmp(simplify,'simplify') 
        for k=1:length(Jdyn_comb(:,1))
            J2dyn(k,:,:) = simplify(J2dyn(k,:,:));
            if ~isempty(J2con)
                J2con(k,:,:) = simplify(J2con(k,:,:));
            end
        end
    elseif strcmp(simplify,'collect')
        for k=1:length(Jdyn_comb(:,1))
            J2dyn(k,:,:) = collect(J2dyn(k,:,:),vars.x);
            if ~isempty(J2con)
                J2con(k,:,:) = collect(J2con(k,:,:),vars.x);
            end
        end
    end
end

% compute third order derivatives
function [J3dyn, J3con] = thirdOrderDerivatives(J2dyn,J2con,vars,simplify)

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

    %compute third order jacobians using 'LR' variables
    % dynamic part
    for k=1:length(J2dyn(:,1,1))
        for l=1:length(J2dyn(1,:,1))
            %Calculate 3rd order Jacobians
            if ~isempty(find(J2dyn(k,l,:), 1))
                J3dyn(k,l,:,:)=jacobian(reshape(J2dyn(k,l,:),[nrOfVars,1]),vars.z);
            end
        end
    end
    %constraint part
    if ~isempty(J2con)
        for k=1:length(J2con(:,1,1))
            for l=1:length(J2con(1,:,1))
                %Calculate 3rd order Jacobians
                if ~isempty(find(J2con(k,l,:), 1))
                    J3con(k,l,:,:)=jacobian(reshape(J2con(k,l,:),[nrOfVars,1]),vars.z);
                end
            end
        end
    end

    % potentially simplify expression            
    if strcmp(simplify,'simplify') 
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
    elseif strcmp(simplify,'collect')
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

%------------- END OF CODE --------------