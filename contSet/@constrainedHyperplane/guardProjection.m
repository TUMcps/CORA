function Rguard = guardProjection(obj,A,B,t_hit,tmin,tmax,R0,options)
% hyperplaneMap - computes the reachable set of the system within a 
% location, detects the guard set that is hit and computes the new set on
% the hyperplane and the subsequent mapping.
%
% Syntax:  
%    [TP,R,activeGuards,Rjump,Rcont] =
%    hyperplaneMap(obj,tStart,R0,options)
%
% Inputs:
%    obj - location object
%    tStart - start time
%    R0 - initial reachable set
%    options - options struct
%
% Outputs:
%    TP - time point struct; e.g. contains time vector of minimum times for reaching guard sets
%    R - cell array of reachable sets
%    activeGuards - active guards
%    Rjump - reachable set after jump according to the reset map
%    Rcont - reachable set due to continuous evolution without guard or
%    invariant inmtersection
%
% Example: 
%
% Other m-files required: initReach, reach, potOut, potInt, guardIntersect,
% reset
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      10-August-2011
% Last update:  30-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

%perform guard projection without considering constraints
Rguard = guardProjection(obj.h,A,B,t_hit,tmin,tmax,R0,options);

%check if constraints are violated
if ~constrSat(Rguard, obj.C, obj.d)
    Pconstr = pplPolytope(obj.C, obj.d);
    Pzono = polytope(Rguard);
    
    %intersect constraint and zonotope
    P = Pconstr & Pzono;
    
    %enclose result by zonotope
    Rguard = zonotope(vertices(P));
end


%------------- END OF CODE --------------