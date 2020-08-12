function [t,x,loc,xJump] = simulate(obj,params)
% simulate - simulates the system within a location, detects the guard set
% that is hit and computes the reset
%
% Syntax:
%    [t,x,loc,xJump] = simulate(obj,params)
%
% Inputs:
%    obj - location object
%    params - struct storing the simulation parameter
%
% Outputs:
%    t - time vector
%    x - state vector
%    loc - next location
%    xJump - state after jump according to the reset map
%
% Example:
%
% Other m-files required: eventFcn, indexList, reset
% Subfunctions: none
% MAT-files required: none
%
% See also: reach

% Author:       Matthias Althoff
% Written:      03-May-2007
% Last update:  10-August-2011
%               13-March-2015
%               17_August-2015
%               10-September-2015
%               19-April-2016
% Last revision:---

%------------- BEGIN CODE --------------

    % get current location
    currentLoc = params.loc;
    params = rmfield(params,'loc');

    % convert all guard sets to halfspace-representation (polytope)
    for i = 1:length(obj.transition)
       obj.transition{i} = guard2polytope(obj.transition{i}); 
    end

    % define event function from halfspace inequalities
    eventOptions = odeset('Events',eventFcn(obj));
    options = odeset(eventOptions);
    
    % simulate continuous dynamics
    [t,x,index] = simulate(obj.contDynamics,params,options);
    
    % determine the guard set which is crossed by the trajectory
    if ~isempty(index)                          % final time not reached
        
        % determine active guard
        [list] = indexList(obj);
        
        nActivatedGuards = length(index);
        activatedGuards = [];
        
        % check if invariant set was left without hitting a guard set
        if all(list(index) == 0)
            error('Trajectory left the invariant set without hitting a guard set!'); 
        end
        
        % loop over all active guards
        for iActivatedGuard = 1:nActivatedGuards
            
            guard = list(index(iActivatedGuard));
            
            if guard ~= 0 && (~any(activatedGuards == guard))
                
                guardSet = obj.transition{guard}.guard;
                
                % check whether only one halfspace has beed crossed or if 
                % the point is indeed inside the guard set
                if in(guardSet, x(end,:)', 1e-6)
                    
                    % next location and reset function
                    loc = obj.transition{guard}.target;
                    xJump=reset(obj.transition{guard},x(end,:));
                    break;
                    
                else
                    
                    % only one halfspace crossed -> continue simulation in
                    % the current location
                    loc = currentLoc;
                    xJump = x(end,:);
                end
                
                activatedGuards = [activatedGuards,guard];
            end
        end
        
    else                                        % final time reached
        loc=[]; xJump=[];
    end

%------------- END OF CODE --------------
