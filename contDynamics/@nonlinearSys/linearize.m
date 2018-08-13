function [obj,linSys,linOptions] = linearize(obj,options,R)
% linearize - linearizes the nonlinear system; linearization error is not
% included yet
%
% Syntax:  
%    [obj,linSys,linOptions] = linearize(obj,options,R)
%
% Inputs:
%    obj - nonlinear system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    obj - linear system object
%    linSys - linear system object
%    linOptions - options for the linearized system
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      29-October-2007 
% Last update:  22-January-2008
%               29-June-2009
%               04-August-2016
%               15-August-2016
%               12-September-2017
% Last revision:---

%------------- BEGIN CODE --------------

%linearization point p.u of the input is the center of the input u
p.u=center(options.U)+options.uTrans;

%obtain linearization point
if isfield(options,'linearizationPoint')
    p.x = options.linearizationPoint;
else
    %linearization point p.x of the state is the center of the last reachable 
    %set R translated by 0.5*f0*delta_t
    t=0; %time invariant system
    f0prev=obj.mFile(t,center(R),p.u);
    try %if time step not yet created
        p.x=center(R)+f0prev*0.5*options.timeStep;
    catch
        disp('time step not yet created');
        p.x=center(R);
    end
end

%substitute p into the system equation in order to obtain the constant
%input
t=0; %time invariant system
f0=obj.mFile(t,p.x,p.u);

%substitute p into the Jacobian with respect to x and u to obtain the
%system matrix A and the input matrix B
[A,B]=obj.jacobian(p.x,p.u);

%set up otions for linearized system
linOptions=options;

linOptions.uTrans=f0; %B*Ucenter from linOptions.U not added as the system is linearized around center(U)
linOptions.U=B*(options.U+(-center(options.U)));
linOptions.originContained=0;

%set up linearized system
linSys = linearSys('linSys',A,1); %B=1 as input matrix encountered in uncertain inputs

%save constant input
obj.linError.f0=f0;

%save linearization point
obj.linError.p=p;


%------------- END OF CODE --------------