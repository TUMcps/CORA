function [R,Rjump_,res] = reach(loc,params,options)
% reach - computes the reachable set of the system within a location and
%    determines the intersections with the guard sets
%
% Syntax:
%    [R,Rjump_,res] = reach(loc,params,options)
%
% Inputs:
%    loc - location object
%    params - model parameters
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
tStart_interval = params.tStart;
params.tStart = infimum(params.tStart);

% adapt specifications
[specReach,specCheck] = aux_adaptSpecs(loc,options.specification);

% since we require the reachable set for the guard intersection and not the
% output set, we set the internal option 'compOutputSet' to false; the
% output set will then be computed in hybridAutomaton/reach after all
% computation in the location are finished
options.compOutputSet = false;

% compute reachable set for the continuous dynamics until the reachable
% set is fully located outside the invariant set
R = reach(loc.contDynamics,params,options,specReach);

% loop over all reachable sets (the number of reachable sets may
% increase if the sets are split during the computation)
for i=1:size(R,1)

    % determine all guard sets of the current location which any
    % reachable set intersects
    [guards,setIndices,setType] = potInt(loc,R(i),params.finalLoc);

    % compute intersections with the guard sets
    [Rguard,actGuards,minInd,maxInd] = ...
            guardIntersect(loc,guards,setIndices,setType,R(i),params,options);

    % compute reset and get target location
    Rjump_ = struct('set',cell(1,0),'time',cell(1,0),...
        'loc',cell(1,0),'parent',cell(1,0));

    for j=1:length(Rguard)
        
        iGuard = actGuards(j);
        
        % compute reset
        Rjump_(j,1).set = evaluate(loc.transition(iGuard).reset,Rguard{j},params.U);
        
        % target location and parent reachable set
        Rjump_(j,1).loc = loc.transition(iGuard).target;
        Rjump_(j,1).parent = R(i).parent + 1;
        
        % time interval for the guard intersection
        if strcmp(setType,'time-interval')
            tMin = infimum(R.timeInterval.time{minInd(j)});
            tMax = supremum(R.timeInterval.time{maxInd(j)}) + 2*rad(tStart_interval);
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
    R(i) = updateTime(R(i),tStart_interval);
    
    % check if specifications are violated
    if ~isempty(specCheck)
        res = check(specCheck,R(i).timeInterval.set{end});
        if ~res
            return; 
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [specReach,specCheck] = aux_adaptSpecs(loc,specs)

% add the invariant as a specification so that reachability analysis of
% continuous dynamics can exit once the reachable set fully leaves the
% invariant
specReach = specification(loc.invariant,'invariant');

% add other specification to the list of specifications: note that unsafe
% sets must intersect the invariant to be valid
if ~isempty(specs)
    % double check because numel(specification()) == 1
    for i=1:numel(specs)
        if ~strcmp(specs(i).type,'unsafeSet') ...
                || isIntersecting(loc.invariant,specs(i).set)
            specReach = add(specReach,specs(i));
        end
    end
end

% for the check, we skip the 'invariant' specification
specCheck = specReach(2:end);

end

% ------------------------------ END OF CODE ------------------------------
