function [t,z,ind] = simulate(obj,params,varargin)
% simulate - simulates the system within a location
%
% Syntax:  
%    [t,z] = simulate(obj,params)
%    [t,z,ind] = simulate(obj,params,options)
%
% Inputs:
%    obj - nonlinDASys object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .x0: initial point
%       .y0guess: 
%       .u: piecewise constant input signal u(t) specified as a matrix
%           for which the number of rows is identical to the number of
%           system input
%    options - ODE45 options (for hybrid systems)
%
% Outputs:
%    t - time vector
%    z - state vector (dynamic and algebraic state variables)
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
% Written:      03-May-2007 
% Last update:  12-March-2008
%               19-August-2016
%               08-May-2020 (MW, update interface)
% Last revision:---

%------------- BEGIN CODE --------------

% specify mass matrix
M = diag([ones(1,obj.dim),zeros(1,obj.nrOfConstraints)]);

% add mass matrix to the options struct
if nargin > 2
   options = varargin{1}; 
   options = odeset(options, 'Mass', M, 'MStateDependence', 'none');
else
   options = odeset('Mass', M, 'MStateDependence', 'none');
end
options = odeset(options,'RelTol',1e-7,'AbsTol',1e-10,'NormControl','on');

% initial state is combination of dynamic and algebraic state variables
if length(params.x0) == (obj.dim + obj.nrOfConstraints)
    z0 = params.x0;
else
    %extract dynamic and algebraic initial state, as well as the input
    y0 = params.y0guess;
    %ensure consistent initial state
    y0 = consistentInitialState(obj, params.x0, y0, params.u);
    %update combined initial state
    z0 = [params.x0;y0];
end

try
    [t,z,~,~,ind] = ode15s(getfcn(obj,params),...
        [params.tStart,params.tFinal],z0,options);
    %[t,z,te,xe,index] = ode23tb(getfcn(obj,params),...
    %   [params.tStart,params.tFinal],z0,options);
catch
    [t,z] = ode15s(getfcn(obj,params),[params.tStart,params.tFinal],z0,options);
    %[t,z] = ode23tb(getfcn(obj,params),[params.tStart,params.tFinal],z0,options);
    ind=[];
end

%------------- END OF CODE --------------