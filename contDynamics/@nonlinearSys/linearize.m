function [obj,linSys,linOptions] = linearize(obj,options,R)
% linearize - linearizes the nonlinear system; linearization error is not
%    included yet
%
% Syntax:
%    [obj,linSys,linOptions] = linearize(obj,options,R)
%
% Inputs:
%    obj - nonlinearSys object
%    options - options struct
%    R - reachable set
%
% Outputs:
%    obj - nonlinearSys object
%    linSys - linearSys object
%    linOptions - options for the linearized system
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       29-October-2007 
% Last update:   22-January-2008
%                29-June-2009
%                04-August-2016
%                15-August-2016
%                12-September-2017
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%linearization point p.u of the input is the center of the input set
p.u = options.uTrans;

%obtain linearization point
if isfield(options,'linearizationPoint')
    p.x = options.linearizationPoint;
elseif isfield(options,'refPoints')
    currentStep = round((options.t-options.tStart)/options.timeStep)+1;
    p.x = 1/2*sum(options.refPoints(:,currentStep:currentStep+1),2);
else
    %linearization point p.x of the state is the center of the last
    %reachable set R translated by 0.5*delta_t*f0
    f0prev=obj.mFile(center(R),p.u);
    try %if time step not yet created
        p.x=center(R)+f0prev*0.5*options.timeStep;
    catch
        disp('time step not yet created');
        p.x=center(R);
    end
end

%substitute p into the system equation to obtain the constant input
f0=obj.mFile(p.x,p.u);

%substitute p into the Jacobian with respect to x and u to obtain the
%system matrix A and the input matrix B
[A,B]=obj.jacobian(p.x,p.u);
A_lin = A;
B_lin = B;
linOptions=options;
if strcmp(options.alg,'linRem')
    %in order to compute dA,dB, we use the reachability set computed
    %for one step in initReach
    [dA,dB] = lin_error2dAB(options.Ronestep,options.U,obj.hessian,p);
    A = matZonotope(A,{dA});
    B = matZonotope(B,{dB});
    linSys = linParamSys(A,1,'constParam');
    linOptions.compTimePoint = 1;
else
    %set up linearized system
    linSys = linearSys('linSys',A,1); %B=1 as input matrix encountered in uncertain inputs
end

%set up options for linearized system
linOptions.U = B*(options.U+options.uTrans-p.u);
Ucenter = center(linOptions.U);
linOptions.U = linOptions.U - Ucenter;
linOptions.uTrans = zonotope([f0 + Ucenter,zeros(size(f0,1),1)]);
%linOptions.uTrans = zonotope(f0 + Ucenter);
linOptions.originContained = false;

%save constant input
obj.linError.f0=f0;

%save linearization point
obj.linError.p=p;

% ------------------------------ END OF CODE ------------------------------
