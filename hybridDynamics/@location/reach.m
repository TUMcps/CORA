function [R,Rjump_,res] = reach(loc,R0,tStart,options)
% reach - computes the reachable set of the system within a location and
%    determines the intersections with the guard sets
%
% Syntax:
%    [R,Rjump_,res] = reach(loc,R0,tStart,options)
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

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       07-May-2007 
% Last update:   17-August-2007
%                31-July-2016
%                19-August-2016
%                09-December-2019 (NK, integrated singleSetReach)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
Rjump = struct('set',cell(1,0),'time',cell(1,0),...
    'loc',cell(1,0),'parent',cell(1,0));

% split options into params and options for contDynamics/reach
[params,options_] = splitIntoParamsOptions(options);
% additional adaptations
params.tStart = infimum(tStart);
params.R0 = R0;

% adapt specifications
spec = specification(loc.invariant,'invariant');

if ~isempty(options.specification)
    spec = add(options.specification,spec);
end

% since we require the reachable set for the guard intersection and not the
% output set, we set the internal option 'compOutputSet' to false; the
% output set will then be computed in hybridAutomaton/reach after all
% computation in the location are finished
options_.compOutputSet = false;

% compute reachable set for the continuous dynamics until the reachable
% set is fully located outside the invariant set
R = reach(loc.contDynamics,params,options_,spec);

% loop over all reachable sets (the number of reachable sets may
% increase if the sets are split during the computation)
for i=1:size(R,1)

    % determine all guard sets of the current location which any
    % reachable set intersects
    [guards,setIndices,setType] = potInt(loc,R(i),options);

    % compute intersections with the guard sets
    [Rguard,actGuards,minInd,maxInd] = ...
            guardIntersect(loc,guards,setIndices,setType,R(i),options);

    % compute reset and get target location
    Rjump_ = struct('set',cell(1,0),'time',cell(1,0),...
        'loc',cell(1,0),'parent',cell(1,0));

    for j=1:length(Rguard)
        
        iGuard = actGuards(j);
        
        % compute reset
        Rjump_(j,1).set = reset(loc.transition(iGuard),Rguard{j},options.U);  
        
        % target location and parent reachable set
        Rjump_(j,1).loc = loc.transition(iGuard).target;
        Rjump_(j,1).parent = R(i).parent + 1;
        
        % time interval for the guard intersection
        if strcmp(setType,'time-interval')
            tMin = infimum(R.timeInterval.time{minInd(j)});
            tMax = supremum(R.timeInterval.time{maxInd(j)}) + 2*rad(tStart);
        else
            tMin = R.timePoint.time{minInd(j)};
            tMax = R.timePoint.time{maxInd(j)};
        end
        
        Rjump_(j,1).time = interval(tMin,tMax);
    end
    
    Rjump = [Rjump;Rjump_];

    % remove the parts of the reachable sets outside the invariant
    if options.intersectInvariant
        R(i) = potOut(loc,R(i),minInd,maxInd);
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

% ------------------------------ END OF CODE ------------------------------
