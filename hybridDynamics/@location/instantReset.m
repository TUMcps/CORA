function [R,Rjump,res] = instantReset(loc,R0,tStart,options)
% instantReset - executes an instant transition, i.e., a transition
%    with a full-dimensional set, leading to a reset at time 0 using the
%    initial set
%
% Syntax:
%    [R,Rjump,res] = instantReset(loc,R0,tStart,options)
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
%    res - true (consistent output argument structure with location/reach)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/reach

% Authors:       Mark Wetzlinger
% Written:       18-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create 'reachable set' from instant reset
timePoint.set = {R0};
timePoint.time = {tStart};
timeInt = [];

% instantiate reachSet object
R = reachSet(timePoint,timeInt);

% reset initial set instantly according to reset function
Rjump{1,1}.set = reset(loc.transition(options.instantTransition),...
                R0,options.U);

% target location of instant transition
Rjump{1,1}.loc = loc.transition(options.instantTransition).target;

% branch number always 1
Rjump{1,1}.parent = 1;

% no time passes
Rjump{1,1}.time = tStart;

% no violation of any specification, since no time has passed
res = true;

% ------------------------------ END OF CODE ------------------------------
