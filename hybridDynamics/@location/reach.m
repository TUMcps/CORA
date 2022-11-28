function [R,Rjump_,res] = reach(loc,R0,tStart,options)
% reach - computes the reachable set of the system within a location and
%    determines the intersections with the guard sets
%
% Syntax:  
%    [R,Rjump,res] = reach(loc,R0,tStart,options)
%
% Inputs:
%    loc - location object
%    R0 - initial reachable set
%    tStart - start time
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - reachable set due to continuous evolution
%    Rjump_ - list of guard set intersections with the corresponding sets
%    res - true/false whether specifications are satisfied
%
% See also: hybridAutomaton/reach

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      07-May-2007 
% Last update:  17-August-2007
%               31-July-2016
%               19-August-2016
%               09-December-2019 (NK, integrated singleSetReach)
% Last revision:---

%------------- BEGIN CODE --------------

res = true;
Rjump = {};

% adapt options
[params,options_] = adaptOptions(loc,options);

params.tStart = infimum(tStart);
params.R0 = R0;

% adapt specifications
spec = specification(loc.invariant,'invariant');

if ~isempty(options.specification)
    spec = add(options.specification,spec);
end

% compute reachable set for the continuous dynamics until the reachable
% set is fully located outside the invariant set
R = reach(loc.contDynamics,params,options_,spec);

% loop over all reachable sets (the number of reachable sets may
% increase if the sets are split during the computation)
for i=1:size(R,1)

    % determine all guard sets of the current location which any
    % reachable set intersects
    [guards,setIndices] = potInt(loc,R(i).timeInterval.set,options);

    % compute intersections with the guard sets
    [Rguard,actGuards,minInd,maxInd] = ...
                    guardIntersect(loc,guards,setIndices,R(i),options);

    % compute reset and get target location
    Rjump_ = cell(length(Rguard),1);

    for j=1:length(Rguard)
        
        iGuard = actGuards(j);
        
        % compute reset
        Rjump_{j,1}.set = reset(loc.transition{iGuard},Rguard{j},options.U);  
        
        % target location and parent reachable set
        Rjump_{j,1}.loc = loc.transition{iGuard}.target;
        Rjump_{j,1}.parent = R(i).parent + 1;
        
        % time interval for the guard intersection
        tMin = infimum(R.timeInterval.time{minInd(j)});
        tMax = supremum(R.timeInterval.time{maxInd(j)}) + 2*rad(tStart);
        
        Rjump_{j,1}.time = interval(tMin,tMax);
    end
    
    Rjump = [Rjump;Rjump_];

    % remove the parts of the reachable sets outside the invariant
    if isfield(options,'intersectInvariant') && options.intersectInvariant
        R(i) = potOut(loc,R(i),minInd,maxInd,options);
    end
    
    % update times of the reachable set due to uncertain initial time
    R(i) = updateTime(R(i),tStart);
    
    % check if specifications are violated
    if ~isempty(options.specification)
        res = check(options.specification,R(i).timeInterval.set{end}); 
        if ~res
            return; 
        end
    end
end

%------------- END OF CODE --------------