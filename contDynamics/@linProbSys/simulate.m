function [obj,t,x,index] = simulate(obj,opt,tstart,tfinal,x0,options)
% simulate - simulates the linear interval system within a location
%
% Syntax:  
%    [t,x,index] = simulate(obj,tstart,tfinal,x0,options)
%
% Inputs:
%    obj - linIntSys object
%    tstart - start time
%    tfinal - final time
%    x0 - initial state 
%    options - contains, e.g. the events when a guard is hit
%
% Outputs:
%    obj - linIntSys object
%    t - time vector
%    x - state vector
%    index - returns the event which has been detected
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      16-May-2007 
% Last update:  26-February-2008
% Last revision:---

%------------- BEGIN CODE --------------

%self-programmed euler solver
h=opt.timeStep/5;
t(1)=tstart;
x(:,1)=x0;

%obtain dimension
dim=length(obj.A);

for i=1:(ceil((tfinal-tstart)/h)+1) %+1 to enforce that simulation is quit
    
    %compute random value from the noise signal
    mu=zeros(dim,1);
    Sigma=1/h*eye(dim);
    u=obj.C*mvnrnd(mu,Sigma)';
    
    %next state
    x(:,i+1)=expm(obj.A*h)*x(:,i)+... %initial solution
        inv(obj.A)*(expm(obj.A*h)-eye(length(obj.A)))*(opt.u+u); %input solution
end
index=[];
x=x';
%ode4: fixed step size Runge-Kutta

    
%------------- END OF CODE --------------