function [nlnsys,linsys,linParams,linOptions] = linearize(nlnsys,R,params,options)
% linearize - linearizes the nonlinear system; linearization error is not
%    included yet
%
% Syntax:
%    [nlnsys,linsys,linParams,linOptions] = linearize(nlnsys,R,params,options)
%
% Inputs:
%    nlnsys - nonlinearSys object
%    R - reachable set (only required if no linearization point given)
%    params - model parameters
%    options - options struct
%
% Outputs:
%    nlnsys - nonlinearSys object
%    linsys - linearSys object, linParamSys object
%    linParams - model parameter for the linearized system
%    linOptions - options for the linearized system
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       29-October-2007 
% Last update:   22-January-2008
%                29-June-2009
%                04-August-2016
%                15-August-2016
%                12-September-2017
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% linearization point p.u of the input is the center of the input set
p.u = params.uTrans;

% obtain linearization point
if isfield(options,'linearizationPoint')
    p.x = options.linearizationPoint;
elseif isfield(options,'refPoints')
    currentStep = round((options.t-params.tStart)/options.timeStep)+1;
    p.x = 1/2*sum(options.refPoints(:,currentStep:currentStep+1),2);
else
    % center of start set
    cx = center(R);
    % linearization point p.x of the state is the center of the last
    % reachable set R translated by 0.5*delta_t*f0
    f0prev = nlnsys.mFile(cx,p.u);
    try %if time step not yet created
        p.x = cx + f0prev*0.5*options.timeStep;
    catch
        disp('time step not yet created');
        p.x = cx;
    end
end

%substitute p into the system equation to obtain the constant input
f0 = nlnsys.mFile(p.x,p.u);

%substitute p into the Jacobian with respect to x and u to obtain the
%system matrix A and the input matrix B
[A,B] = nlnsys.jacobian(p.x,p.u);
linOptions=options;
if strcmp(options.alg,'linRem')
    %in order to compute dA,dB, we use the reachability set computed
    %for one step in initReach
    [dA,dB] = lin_error2dAB(options.Ronestep,params.U,nlnsys.hessian,p);
    A = matZonotope(A,dA);
    B = matZonotope(B,dB);
    linsys = linParamSys(A,1,'constParam');
    linOptions.compTimePoint = true;
else
    % set up linearized system (B = 1 as input matrix encountered in
    % uncertain inputs)
    linsys = linearSys([nlnsys.name '_linearized'],A,1);
end

%set up options for linearized system
linParams.U = B*(params.U+params.uTrans-p.u);
Ucenter = center(linParams.U);
linParams.U = linParams.U - Ucenter;
linParams.uTrans = zonotope(f0 + Ucenter,zeros(size(f0,1),1));
%linParams.uTrans = zonotope(f0 + Ucenter);
linOptions.originContained = false;

%save constant input
nlnsys.linError.f0=f0;

%save linearization point
nlnsys.linError.p=p;

% ------------------------------ END OF CODE ------------------------------
