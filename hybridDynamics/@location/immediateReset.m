function [R,Rjump,res] = immediateReset(loc,R0,tStart,options)
% immediateReset - executes an immediate transition, i.e., a transition
%    with no guard set where the initial set is reset at time 0
%
% Syntax:  
%    [R,Rjump,res] = immediateReset(loc,R0,tStart,options)
%
% Inputs:
%    loc - location object
%    R0 - initial reachable set
%    tStart - start time
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - reachable set due to continuous evolution
%    Rjump - list of guard set intersections with the corresponding sets
%    res - false (consistent output argument structure with location/reach)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/reach

% Author:       Mark Wetzlinger
% Written:      18-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create 'reachable set' from immediate reset
timePoint.set = {R0};
timePoint.time = {tStart};
timeInt.set = {};
timeInt.time = {};

% instantiate reachSet object
R = reachSet(timePoint,timeInt);

% reset initial set immediately according to reset function
Rjump{1,1}.set = reset(loc.transition{options.immediateTransition},...
                R0,options.U);

% target location of immediate transition
Rjump{1,1}.loc = loc.transition{options.immediateTransition}.target;

% branch number always 1
Rjump{1,1}.parent = 1;

% no time passes
Rjump{1,1}.time = tStart;

% no violation of any specification, since no time has passed
res = true;

%------------- END OF CODE --------------