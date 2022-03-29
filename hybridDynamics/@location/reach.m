function [R,Rjump_,res] = reach(obj,R0,tStart,options)
% reach - computes the reachable set of the system within a location and
%         determines the intersections with the guard sets
%
% Syntax:  
%    [R,Rjump] = reach(obj,R0,tStart,options)
%
% Inputs:
%    obj - location object
%    R0 - initial reachable set
%    tStart - start time
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - reachable set due to continuous evolution
%    Rjump - list of guard set intersections with the corresponding sets
%
% See also: hybridAutomaton/reach

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      07-May-2007 
% Last update:  17-August-2007
%               31-July-2016
%               19-August-2016
%               09-December-2019 (NK, integrated singleSetReach)
% Last revision: ---

%------------- BEGIN CODE --------------

    res = 1;
    Rjump = {};

    % adapt options
    [params,options_] = adaptOptions(obj,options);
    
    params.tStart = infimum(tStart);
    params.R0 = R0;
    
    % adapt specifications
    spec = specification(obj.invariant,'invariant');
    
    if ~isempty(options.specification)
        spec = add(options.specification,spec);
    end
    
    % compute reachable set for the continuous dynamics until the reachable
    % set is fully located outside the invariant set
    R = reach(obj.contDynamics,params,options_,spec);
    
    % loop over all reachable sets
    for i = 1:size(R,1)

        % determine all guard sets which the reachable set intersects
        [guards,setIndices] = potInt(obj,R(i).timeInterval.set,options);

        % compute intersections with the guard sets
        [Rguard,actGuards,minInd,maxInd] = ....
                        guardIntersect(obj,guards,setIndices,R(i),options);

        % compute reset and get target location
        Rjump_ = cell(length(Rguard),1);

        for j = 1:length(Rguard)
             
            iGuard = actGuards(j);
            
            % compute reset
            Rjump_{j,1}.set = reset(obj.transition{iGuard},Rguard{j},options.U);  
            
            % target location and parent reachable set
            Rjump_{j,1}.loc = obj.transition{iGuard}.target;
            Rjump_{j,1}.parent = R(i).parent + 1;
            
            % time interval for the guard intersection
            tMin = infimum(R.timeInterval.time{minInd(j)});
            tMax = supremum(R.timeInterval.time{maxInd(j)}) + 2*rad(tStart);
            
            Rjump_{j,1}.time = interval(tMin,tMax);
        end
        
        Rjump = [Rjump;Rjump_];

        % remove the parts of the reachable sets outside the invariant
        if isfield(options,'intersectInvariant') && options.intersectInvariant
            R(i) = potOut(obj,R(i),minInd,maxInd,options);
        end
        
        % update times of the reach. set due to uncertain initial time
        R(i) = updateTime(R(i),tStart);
        
        % check if specifications are violated
        if ~isempty(options.specification)
           res = check(options.specification,R(i).timeInterval.set{end}); 
           if ~res
              return; 
           end
        end
    end
end


% Auxiliary Functions -----------------------------------------------------

function R = updateTime(R,time)
% update the times of the reachable set due to the uncertain initial time

    deltaT = time - infimum(time);
    timePoint = R.timePoint;
    timeInt = R.timeInterval;
    
    for i = 1:length(timePoint.time)
       timePoint.time{i} = timePoint.time{i} + deltaT; 
       timeInt.time{i} = timeInt.time{i} + deltaT; 
    end
    
    R = reachSet(timePoint,timeInt,R.parent);
end

%------------- END OF CODE --------------