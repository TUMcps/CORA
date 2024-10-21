function [obj,linsys,linParams,linOptions] = linearize(obj,R,params,options)
% linearize - linearizes the nonlinear system with uncertain parameters;
%    linearization error is not included yet
%
% Syntax:
%    [obj,linsys,linParams,linOptions] = linearize(obj,R,params,options)
%
% Inputs:
%    obj - nonlinParamSys object
%    R - actual reachable set
%    params - model parameters
%    options - options struct
%
% Outputs:
%    obj - nonlinParamSys object
%    linsys - linParamSys/linearSys object
%    linOptions - options for the linearized system
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-January-2008 
% Last update:   30-June-2009
%                27-May-2011
%                07-June-2017
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


%linearization point p.u of the input is the center of the input u
p.u=center(params.U) + params.uTrans;

% one paramInt for each reachStep
if isfield(options,'paramInts')
    currentStep = round((options.t-params.tStart)/options.timeStep)+1;
    params.paramInt = options.paramInts(:,currentStep);
end

%linearization point p.p of the parameters is the center of the parameter
%set paramInt;
if isa(params.paramInt,'interval')
    p.p = center(params.paramInt);
else
    p.p = params.paramInt; 
end

%linearization point p.x of the state is the center of the last reachable 
%set R translated by f0*delta_t
if isfield(options,'refPoints')
    currentStep = round((options.t-params.tStart)/options.timeStep)+1;
    p.x = 1/2*sum(options.refPoints(:,currentStep:currentStep+1),2);
else
    f0prev=obj.mFile(center(R),p.u,p.p);
    p.x=center(R)+f0prev*0.5*options.timeStep;
end

%two options: uncertain parameters with fixed bounds or fixed parameters
%that can change during the reachability analysis
%option 1: uncertain paramters
if isa(params.paramInt,'interval')

    % obtain cell of possible derivatives
    f0_3D = obj.parametricDynamicFile(p.x,p.u);
    f0cell = arrayfun(@(i) f0_3D(:,:,i),1:obj.nrOfParam+1,...
        'UniformOutput',false);

    %normalize cells
    for i=2:length(f0cell)
        f0cell{1} = f0cell{1} + center(params.paramInt(i-1))*f0cell{i};
        f0cell{i} = rad(params.paramInt(i-1))*f0cell{i};
    end

    %create constant input zonotope
    f0Mat(length(f0cell{1}(:,1)),length(f0cell)) = 0; %init
    for i=1:length(f0cell)
        f0Mat(:,i) = f0cell{i};
    end
    zf = zonotope(f0Mat);

    % substitute p into the Jacobian with respect to x and u to obtain the
    % system matrix A and the input matrix B
    [A,B] = obj.jacobian(p.x,p.u);

    %create matrix zonotopes
    zA = matZonotope(A(:,:,1),A(:,:,2:end));
    zB = matZonotope(B(:,:,1),B(:,:,2:end));
    
    %check if remainder should be included in system matrices
    %(options.alg = 'linRem')
    linOptions = options;
    if strcmp(options.alg,'linRem')
        %in order to compute dA,dB, we use the reachability set computed
        %for one step in initReach
        [dA,dB] = lin_error2dAB(options.Ronestep,params.U,obj.hessian,p,params.paramInt);
        zA = matZonotope(zA.C,cat(3,zA.G,dA));
        zB = matZonotope(zB.C,cat(3,zB.G,dB));
        zf = zonotope(zf.c, [zf.G,zeros(dim(zf),1)]);%needed for dependentHomoSol
        linOptions.compTimePoint = true;
    end
    %set up otions for linearized system
    U = zB*(params.U + params.uTrans + (-p.u));
    Ucenter = center(U);
    linParams.U = U + (-Ucenter);
    linParams.uTrans = Ucenter;
    linParams.Uconst = zf + Ucenter;
    linOptions.originContained = false;

    % instantiate linearized system
    if ~obj.constParam
        % time-variant linear parametric system
        % B=1 as input matrix encountered in uncertain inputs
        linsys = linParamSys(zA,1,'varParam');
    else
        % time-invariant linear parametric system
        % B=1 as input matrix encountered in uncertain inputs
        linsys = linParamSys(zA,1,'constParam');
    end

    %save constant input
    obj.linError.f0=zf; 
    
%option 2: fixed parameters that can change
elseif isnumeric(params.paramInt)
    % obtain derivative
    f0 = obj.mFile(p.x,p.u,p.p);
    
    %substitute p into the Jacobian with respect to x and u to obtain the
    %system matrix A and the input matrix B
    [A,B]=obj.jacobian_freeParam(p.x,p.u,p.p);

    %set up otions for linearized system
    linOptions=options;

    linParams.uTrans=f0; %B*Ucenter from linOptions.U not added as the system is linearized around center(U)
    linParams.U=B*(params.U+(-center(params.U)));
    linOptions.originContained=false;

    %set up linearized system
    linsys = linearSys('linsys',A,1); %B=1 as input matrix encountered in uncertain inputs

    %save constant input
    obj.linError.f0=f0;
end

%save linearization point
obj.linError.p=p;

% ------------------------------ END OF CODE ------------------------------
