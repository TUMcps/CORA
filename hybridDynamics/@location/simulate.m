function [t,x,nextloc,xJump] = simulate(loc,params)
% simulate - simulates the system within a location, detects the guard set
%    that is hit and computes the reset
%
% Syntax:
%    [t,x,nextloc,xJump] = simulate(loc,params)
%
% Inputs:
%    loc - location object
%    params - struct storing the simulation parameter
%
% Outputs:
%    t - time vector
%    x - state vector
%    nextloc - next location
%    xJump - state after jump according to the reset map
%
% Example:
%    ---
%
% Other m-files required: eventFcn, indexList, reset
% Subfunctions: none
% MAT-files required: none
%
% See also: reach

% Authors:       Matthias Althoff
% Written:       03-May-2007
% Last update:   10-August-2011
%                13-March-2015
%                17-August-2015
%                10-September-2015
%                19-April-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get current location
currentLoc = params.loc;
params = rmfield(params,'loc');

% convert all guard sets to halfspace-representation (polytope)
for i=1:length(loc.transition)
    loc.transition(i) = guard2polytope(loc.transition(i)); 
end

% define event function from halfspace inequalities
eventOptions = odeset('Events',eventFcn(loc));
options = odeset(eventOptions);

% simulate continuous dynamics
[t,x,index] = simulate(loc.contDynamics,params,options);

% determine the guard set which is crossed by the trajectory
if ~isempty(index)
    % final time not reached

    % determine active guard
    list = indexList(loc);
    
    nActivatedGuards = length(index);
    activatedGuards = [];
    
    % check if invariant set was left without hitting a guard set
    if all(list(index) == 0)
        throw(CORAerror('CORA:specialError',...
            'Trajectory left the invariant set without hitting a guard set!')); 
    end
    
    % loop over all active guards
    for iActivatedGuard = 1:nActivatedGuards
        
        guard = list(index(iActivatedGuard));
        
        if guard ~= 0 && (~any(activatedGuards == guard))
            
            guardSet = loc.transition(guard).guard;
            
            % check whether only one halfspace has beed crossed or if 
            % the point is indeed inside the guard set
            if contains_(guardSet,x(end,:)','exact',1e-6)
                
                % next location and reset function
                nextloc = loc.transition(guard).target;
                xJump = reset(loc.transition(guard),x(end,:));
                break;
                
            else
                
                % only one halfspace crossed -> continue simulation in
                % the current location
                nextloc = currentLoc;
                xJump = x(end,:);
            end
            
            activatedGuards = [activatedGuards,guard];
        end
    end
    
else
    % final time reached
    nextloc = []; xJump = [];
end

% ------------------------------ END OF CODE ------------------------------
