function [t,x] = simulate(obj,options)
% simulate - simulates the linear interval system within a location
%
% Syntax:  
%    [t,x,ind] = simulate(obj,params,options)
%
% Inputs:
%    obj - linProbSys object
%    options - struct containing the parameters for the simulation
%       .tStart initial time t0
%       .tFinal final time tf
%       .x0 initial point x0
%       .u piecewise constant input signal u(t) specified as a matrix
%           for which the number of rows is identical to the number of
%           system input
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - returns the event which has been detected
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
%               17-July-2020
% Last revision:---

%------------- BEGIN CODE --------------

%self-programmed euler solver
h = options.timeStep/5;
nrOfTimeSteps = ceil((options.tFinal-options.tStart)/h);
t = linspace(options.tStart, options.tFinal, nrOfTimeSteps);
x(:,1) = options.x0;

%obtain dimension
dim=length(obj.A);

for i=1:nrOfTimeSteps-1
    
    %compute random value from the noise signal
    mu=zeros(dim,1);
    Sigma=1/h*eye(dim);
    u=obj.C*mvnrnd(mu,Sigma)';
    
    %next state
    x(:,i+1)=expm(obj.A*h)*x(:,i)+... %initial solution
        inv(obj.A)*(expm(obj.A*h)-eye(length(obj.A)))*(options.u+u); %input solution
end
x=x';

%------------- END OF CODE --------------