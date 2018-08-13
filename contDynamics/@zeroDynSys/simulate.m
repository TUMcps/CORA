function [obj,t,x,index] = simulate(obj,opt,tstart,tfinal,x0,options)
% simulate - simulates the system within a location
%
% Syntax:  
%    [t,x,index] = simulate(obj,tstart,tfinal,x0,options)
%
% Inputs:
%    obj - linearSys object
%    tstart - start time
%    tfinal - final time
%    x0 - initial state 
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

% Author: Matthias Althoff
% Written: 02-November-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

t=opt.tFinal;
x=x0';
te=[];
xe=[];
index=[];

%------------- END OF CODE --------------