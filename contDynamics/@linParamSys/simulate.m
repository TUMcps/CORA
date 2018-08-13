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

% Author: Matthias Althoff
% Written: 16-May-2007 
% Last update: 07-January-2009
% Last revision: ---

%------------- BEGIN CODE --------------

%sample system matrix
tmp = randomSampling(obj.A,1,opt);
obj.sampleMatrix.A = tmp{1};

if isempty(options.Events)
    [t,x] = ode45(getfcn(obj,opt),[tstart, tfinal],x0,options);
else
    [t,x,te,xe,index] = ode45(getfcn(obj,opt),[tstart, tfinal],x0,options);
end
    
    
%------------- END OF CODE --------------