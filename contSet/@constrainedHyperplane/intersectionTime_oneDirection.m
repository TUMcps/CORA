function [tmin,tmax,t_total,Rlast] = intersectionTime_oneDirection(obj,R0,contDynamics,options)
% intersectionTime_oneDirection - computes the time to intersect the
% halfspace in the given flow direction
%
% Syntax:  
%    [tmin,tmax,t_total,Rlast] = intersectionTime_oneDirection(obj,R0,contDynamics,options)
%
% Inputs:
%    obj - halfspace object
%    R0 - initial set
%    contDynamics - continuous dynamics
%    options - options struct
%
% Outputs:
%    tmin - minimum time
%    tmax - maximum time
%    t_total - time elapsed for complete intersection starting at R0
%    Rlast - last reachable set before first intersection with halfspace
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      10-March-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%apply the same method as for the halfspace
[tmin,tmax,t_total,Rlast] = intersectionTime_oneDirection(obj.h,R0,contDynamics,options);



%------------- END OF CODE --------------
