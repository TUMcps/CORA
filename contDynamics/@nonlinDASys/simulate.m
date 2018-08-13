function [obj,t,z,index] = simulate(obj,opt,tstart,tfinal,x0,options)
% simulate - simulates the system within a location
%
% Syntax:  
%    [t,x,index] = simulate(obj,tstart,tfinal,x0,options)
%
% Inputs:
%    obj - linearSys object
%    tstart - start time
%    tfinal - final time
%    z0 - initial state of dynamic and algebraic state variables
%    options - contains, e.g. the events when a guard is hit
%
% Outputs:
%    obj - linearSys object
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
% Written:      03-May-2007 
% Last update:  12-March-2008
%               19-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%specify mass matrix
M = diag([ones(1,obj.dim),zeros(1,obj.nrOfConstraints)]);

%add mass matrix to the options struct
options = odeset(options, 'Mass', M, 'MStateDependence', 'none');
options = odeset(options,'RelTol',1e-7,'AbsTol',1e-10,'NormControl','on');

if length(x0) == (obj.dim + obj.nrOfConstraints) % initial state is combination of dynamic and algebraic state variables
    z0 = x0;
else
    %extract dynamic and algebraic initial state, as well as the input
    y0 = opt.y0guess;
    %ensure consisten initial state
    y0 = consistentInitialState(obj, x0, y0, opt.u);
    %update combined initial state
    z0 = [x0;y0];
end

try
    [t,z,te,xe,index] = ode15s(getfcn(obj,opt),[tstart, tfinal],z0,options);
    %[t,z,te,xe,index] = ode23tb(getfcn(obj,opt),[tstart, tfinal],z0,options);
catch
    [t,z] = ode15s(getfcn(obj,opt),[tstart, tfinal],z0,options);
    %[t,z] = ode23tb(getfcn(obj,opt),[tstart, tfinal],z0,options);
    index=[];
end

%------------- END OF CODE --------------