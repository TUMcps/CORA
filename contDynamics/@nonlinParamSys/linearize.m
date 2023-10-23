function [obj,linSys,linOptions] = linearize(obj,options,R)
% linearize - linearizes the nonlinear system with uncertain parameters;
%    linearization error is not included yet
%
% Syntax:
%    [obj,linSys,linOptions] = linearize(obj,options,R)
%
% Inputs:
%    obj - nonlinParamSys object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    obj - nonlinParamSys object
%    linSys - linParamSys/linearSys object
%    linOptions - options for the linearized system
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       18-January-2008 
% Last update:   30-June-2009
%                27-May-2011
%                07-June-2017
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


%linearization point p.u of the input is the center of the input u
p.u=center(options.U) + options.uTrans;

% one paramInt for each reachStep
if isfield(options,'paramInts')
    currentStep = round((options.t-options.tStart)/options.timeStep)+1;
    options.paramInt = options.paramInts(:,currentStep);
end

%linearization point p.p of the parameters is the center of the parameter
%set paramInt;
if isa(options.paramInt,'interval')
    p.p = center(options.paramInt);
else
    p.p = options.paramInt; 
end

%linearization point p.x of the state is the center of the last reachable 
%set R translated by f0*delta_t
if isfield(options,'refPoints')
    currentStep = round((options.t-options.tStart)/options.timeStep)+1;
    p.x = 1/2*sum(options.refPoints(:,currentStep:currentStep+1),2);
else
    f0prev=obj.mFile(center(R),p.u,p.p);
    p.x=center(R)+f0prev*0.5*options.timeStep;
end

%two options: uncertain parameters with fixed bounds or fixed parameters
%that can change during the reachability analysis
%option 1: uncertain paramters
if isa(options.paramInt,'interval')

    % obtain cell of possible derivatives
    f0cell=obj.parametricDynamicFile(p.x,p.u);

    %normalize cells
    for i=2:length(f0cell)
        f0cell{1} = f0cell{1} + center(options.paramInt(i-1))*f0cell{i};
        f0cell{i} = rad(options.paramInt(i-1))*f0cell{i};
    end

    %create constant input zonotope
    f0Mat(length(f0cell{1}(:,1)),length(f0cell)) = 0; %init
    for i=1:length(f0cell)
        f0Mat(:,i) = f0cell{i};
    end
    zf = zonotope(f0Mat);

    %substitute p into the Jacobian with respect to x and u to obtain the
    %system matrix A and the input matrix B
    [A,B]=obj.jacobian(p.x,p.u);

    %create matrix zonotopes
    zA = matZonotope(A{1},A(2:end));
    zB = matZonotope(B{1},B(2:end));
    
    %check if remainder should be included in system matrices
    %(options.alg = 'linRem')
    linOptions = options;
    if strcmp(options.alg,'linRem')
        %in order to compute dA,dB, we use the reachability set computed
        %for one step in initReach
        [dA,dB] = lin_error2dAB(options.Ronestep,options.U,obj.hessian,p,options.paramInt);
        Z_A = zA.generator;
        Z_A{end+1} = dA;
        zA = matZonotope(zA.center,Z_A);
        Z_B = zB.generator;
        Z_B{end+1} = dB;
        zB = matZonotope(zB.center,Z_B);
        zf = zonotope(zf.c, [zf.G,zeros(dim(zf),1)]);%needed for dependentHomoSol
        linOptions.compTimePoint = true;
    end
    %set up otions for linearized system
    U = zB*(options.U + options.uTrans + (-p.u));
    Ucenter = center(U);
    linOptions.U = U + (-Ucenter);
    linOptions.uTrans = Ucenter;
    linOptions.Uconst = zf + Ucenter;
    linOptions.originContained = false;

    % instantiate linearized system
    if ~obj.constParam
        % time-variant linear parametric system
        % B=1 as input matrix encountered in uncertain inputs
        linSys = linParamSys(zA,1,'varParam');
    else
        % time-invariant linear parametric system
        % B=1 as input matrix encountered in uncertain inputs
        linSys = linParamSys(zA,1,'constParam');
    end

    %save constant input
    obj.linError.f0=zf; 
    
%option 2: fixed parameters that can change
elseif isnumeric(options.paramInt)
    % obtain derivative
    f0 = obj.mFile(p.x,p.u,p.p);
    
    %substitute p into the Jacobian with respect to x and u to obtain the
    %system matrix A and the input matrix B
    [A,B]=obj.jacobian_freeParam(p.x,p.u,p.p);

    %set up otions for linearized system
    linOptions=options;

    linOptions.uTrans=f0; %B*Ucenter from linOptions.U not added as the system is linearized around center(U)
    linOptions.U=B*(options.U+(-center(options.U)));
    linOptions.originContained=false;

    %set up linearized system
    linSys = linearSys('linSys',A,1); %B=1 as input matrix encountered in uncertain inputs

    %save constant input
    obj.linError.f0=f0;
end


%save linearization point
obj.linError.p=p;


% ------------------------------ END OF CODE ------------------------------
