function [obj] = symbolicDerivation(obj,options)
% symbolicDerivation - computes the derivatives of the nonlinear system 
% in a symbolic way; the result is also stored in a m-file and passed by 
% a handle
%
% Syntax:  
%    [obj] = symbolicDerivation(obj)
%
% Inputs:
%    obj - nonlinear parameter system object
%
% Outputs:
%    obj - nonlinear parameter system object
%
% Example: 
%    Text for example...
%
% Other m-files required: ---
% Subfunctions: ---
% MAT-files required: ---
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      25-May-2011
% Last update:  02-June-2017
%               06-July-2017
%               15-July-2017 (NK)
%               16-July-2017 (NK)
% Last revision:---

%------------- BEGIN CODE --------------

%create path
path = [coraroot '/contDynamics/stateSpaceModels'];

%create symbolic variables
[vars]=symVariables(obj,'LRbrackets');

%insert symbolic variables into the system equations
fdyn = obj.mFile(vars.x,vars.u,vars.p);

% check if files should be symbolically generated or not-------------------
% check if dynamics files already exists
generateFiles = 1;
filepathOld = [path filesep obj.name '_lastVersion.mat'];

if exist(filepathOld,'file')
    % load stored dynamics 
    load(filepathOld);
    
    % compare stored and actual dynamics (given in symbolic variables)
    if exist('storedData','var')
        if isequal(fdyn,storedData.fdyn) && options.tensorOrder == storedData.tensorOrder

            % check if the same value for options.simplify has been used
            if isfield(options,'simplify')
                if strcmp(options.simplify,storedData.simplify)
                    generateFiles = 0;
                end
            else
                if strcmp(storedData.simplify,'none')
                    generateFiles = 0; 
                end
            end

            % check if parameters are constant or variable
            if isfield(storedData,'paramConst')
                if isa(options.paramInt,'interval') && storedData.paramConst == 1           
                    generateFiles = 1;
                end
            else
                generateFiles = 1;
            end
        end
    end
end

%--------------------------------------------------------------------------

% create files
if generateFiles

    % store the actual dynamics (given in symbolic variables) and the value 
    % for options.simplify in the /StateSpaceModel-directory inside cora
    storedData.fdyn = fdyn;
    storedData.tensorOrder = options.tensorOrder;
    if isfield(options,'simplify')
        storedData.simplify = options.simplify;
    else
        storedData.simplify = 'none';
    end
    if isa(options.paramInt,'interval')
        storedData.paramConst = 0;
    else
        storedData.paramConst = 1;
    end

    disp('create Jacobians');

    %compute jacobian with respect to the state
    Jdyn.x = jacobian(fdyn,vars.x);
    if isfield(options,'simplify')              %potentially simplify expression
        if strcmp(options.simplify,'simplify') 
            Jdyn.x = simplify(Jdyn.x);
        elseif strcmp(options.simplify,'collect')
            Jdyn.x = collect(Jdyn.x,vars.x);
        else
            error('Wrong value for opions.simplify!. Only values "simplify" and "collect" are valid!')
        end
    end

    %compute jacobian with respect to the input
    Jdyn.u = jacobian(fdyn,vars.u);
    if isfield(options,'simplify')              %potentially simplify expression
        if strcmp(options.simplify,'simplify') 
            Jdyn.u = simplify(Jdyn.u);
        elseif strcmp(options.simplify,'collect')
            Jdyn.u = collect(Jdyn.u,vars.x);
        else
            error('Wrong value for opions.simplify!. Only values "simplify" and "collect" are valid!')
        end
    end

    %store jacobians with respect to parameters (assumption: expressions are linear in parameters)
    %init
    Jp.x=cell(1,obj.nrOfParam+1);
    Jp.u=cell(1,obj.nrOfParam+1);


    %part without parameters
    try
        Jp.x{1} = subs(Jdyn.x,p,zeros(obj.nrOfParam,1));
        Jp.u{1} = subs(Jdyn.u,p,zeros(obj.nrOfParam,1));
        %part with parameters
        I = eye(obj.nrOfParam); %identity matrix
        for i=1:obj.nrOfParam
            Jp.x{i+1} = subs(Jdyn.x,p,I(:,i)) - Jp.x{1};
            Jp.u{i+1} = subs(Jdyn.u,p,I(:,i)) - Jp.u{1};
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

    disp('create 2nd order Jacobians');

    %compute second order jacobians using 'LR' variables
    Jdyn_comb = jacobian(fdyn,[vars.x;vars.u]);
    %dynamic jacobians
    for k=1:length(Jdyn_comb(:,1))
        %Calculate 2nd order Jacobians
        J2dyn(k,:,:)=jacobian(Jdyn_comb(k,:),[vars.x;vars.u]);
    end
    %potentially simplify expression
    if isfield(options,'simplify')             
        if strcmp(options.simplify,'simplify') 
            for k=1:length(Jdyn_comb(:,1))
                J2dyn(k,:,:) = simplify(J2dyn(k,:,:));
            end
        elseif strcmp(options.simplify,'collect')
            for k=1:length(Jdyn_comb(:,1))
                J2dyn(k,:,:) = collect(J2dyn(k,:,:),x);
            end
        else
            error('Wrong value for opions.simplify!. Only values "simplify" and "collect" are valid!')
        end
    end

    %group expressions by parameters
    %init
    Jpxu=cell(1,obj.nrOfParam+1);
    JxuAll=J2dyn;
    %part without parameters
    Jpxu{1} = subs(JxuAll,vars.p,zeros(obj.nrOfParam,1));
    %part with parameters
    I = eye(obj.nrOfParam); %identity matrix
    for i=1:obj.nrOfParam
        Jpxu{i+1} = subs(JxuAll,vars.p,I(:,i)) - Jpxu{1};
    end
    obj.derivative.secondOrderParam=Jpxu;


   if options.tensorOrder>2

        disp('create 3rd order Jacobians');
        dim = length(J2dyn(:,1,1));
        nrOfVars = length(J2dyn(1,:,1));
        J3dyn = sym(zeros(dim,nrOfVars,nrOfVars,nrOfVars));

        %compute third order jacobians using 'LR' variables
        for k=1:length(J2dyn(:,1,1))
            for l=1:length(J2dyn(1,:,1))
                %Calculate 3rd order Jacobians
                if ~isempty(find(J2dyn(k,l,:), 1))
                    J3dyn(k,l,:,:)=jacobian(reshape(J2dyn(k,l,:),[nrOfVars,1]),[vars.x;vars.u]);
                end
            end
        end
        %potentially simplify expression
        if isfield(options,'simplify')              
            if strcmp(options.simplify,'simplify') 
                for k=1:length(J2dyn(:,1,1))
                    for l=1:length(J2dyn(1,:,1))
                        if ~isempty(find(J2dyn(k,l,:), 1))
                            J3dyn(k,l,:,:) = simplify(J3dyn(k,l,:,:));
                        end
                    end
                end
            elseif strcmp(options.simplify,'collect')
                for k=1:length(J2dyn(:,1,1))
                    for l=1:length(J2dyn(1,:,1))
                        if ~isempty(find(J2dyn(k,l,:), 1))
                            J3dyn(k,l,:,:) = collect(J3dyn(k,l,:,:),x);
                        end
                    end
                end
            else
                error('Wrong value for opions.simplify!. Only values "simplify" and "collect" are valid!')
            end
        end
    end

    disp('create mFiles')

    %generate mFile that computes the lagrange remainder and the jacobian
    if ~isempty(Jp)
        createJacobianFile(Jp,path,obj.name);
        createParametricDynamicFile(obj);
    end
    createJacobianFile_freeParam(Jdyn,path,obj.name);
    createRemainderFile(obj,Jpxu,path,obj.name);

    if options.tensorOrder>2
        createHessianTensorFile(J2dyn,path,obj.name,0);
        create3rdOrderTensorFile(J3dyn,path,obj.name);
    else
        createHessianTensorFile(J2dyn,path,obj.name,1);
    end

    % rehash the folder so that new generated files are used
    rehash path;

    % save data so that symbolic computations do not have to be re-computed
    save(filepathOld,'storedData');
end

%------------- END OF CODE --------------