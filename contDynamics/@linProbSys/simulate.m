function [t,x,ind,y] = simulate(obj,params,options)
% simulate - simulates the linear interval system within a location
%
% Syntax:
%    [t,x] = simulate(obj,params,options)
%    [t,x,ind] = simulate(obj,params,options)
%    [t,x,ind,y] = simulate(obj,params,options)
%
% Inputs:
%    obj - linProbSys object
%    params - struct containing the parameters for the simulation
%       .tStart initial time t0
%       .tFinal final time tf
%       .x0 initial point x0
%       .u piecewise constant input signal u(t) specified as a matrix
%           for which the number of rows is identical to the number of
%           system input
%    options - settings
%       .timeStep time step size
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - returns the event which has been detected
%    y - output vector
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       16-May-2007 
% Last update:   26-February-2008
%                17-July-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargout == 3
    ind = [];
end
if nargout == 4
    CORAwarning('CORA:contDynamics',"Output trajectories not supported for class linProbSys!");
    y = [];
end

%self-programmed euler solver
h = options.timeStep/5;
nrOfTimeSteps = ceil((params.tFinal-params.tStart)/h);
t = linspace(params.tStart, params.tFinal, nrOfTimeSteps);
x(:,1) = params.x0;

%obtain dimension
n = obj.nrOfDims;

for i=1:nrOfTimeSteps-1
    % compute random value from the noise signal
    mu = zeros(n,1);
    Sigma = 1/h*eye(n);
    u = obj.C*mvnrnd(mu,Sigma)';
    
    % next state: initial solution + input solution
    x(:,i+1) = expm(obj.A*h)*x(:,i) + ...
        inv(obj.A)*(expm(obj.A*h) - eye(length(obj.A)))*(params.u+u);
end
x=x';

% ------------------------------ END OF CODE ------------------------------
