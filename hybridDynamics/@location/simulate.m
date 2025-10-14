function [t,x,nextloc,xJump] = simulate(loc,params)
% simulate - simulates the system within a location, detects the guard set
%    that is hit and computes the reset
%
% Syntax:
%    [t,x,nextloc,xJump] = simulate(loc,params)
%
% Inputs:
%    loc - location object
%    params - struct storing the simulation parameters
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
% Other m-files required: eventFcn, priv_indexList, reset
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

if size(params.x0,2) > size(params.x0,1)
    params.x0 = params.x0';
end

% check if initial state is inside the invariant of the current location
if ~contains_(loc.invariant,params.x0,'exact',1e-6,0,false,false)
    throw(CORAerror('CORA:specialError',...
            'Trajectory is located outside the invariant set after transition!')); 
end

% convert all guard sets to halfspace-representation (polytope)
for i=1:length(loc.transition)
    loc.transition(i) = guard2polytope(loc.transition(i)); 
end

% define event function from halfspace inequalities
fun = eventFcn(loc);
eventOptions = odeset('Events',fun);
options = odeset(eventOptions);

% simulate continuous dynamics
[t,x,index] = simulate(loc.contDynamics,params,options);

% determine the guard set which is crossed by the trajectory
if isempty(index)
    % final time reached
    nextloc = []; xJump = [];
    return
end

% final time not reached

% re-compute index because simulate often just returns the first of 
% potentially multiple events that occurred at the same time
indexNew = find(abs(fun(t(end),x(:,end))) < 1e-6);
index = unique([index;indexNew]);

% determine active guard
list = priv_indexList(loc);

nActivatedGuards = length(index);
activatedGuards = [];

% check if invariant set was left without hitting a guard set
if all(list(index) == 0)
    throw(CORAerror('CORA:specialError',...
        'Trajectory left the invariant set without hitting a guard set!')); 
end

% do not consider events where the state jumps from outside of the
% invariant back inside (it can happen due to numeric imprecisions that the
% state is slightly located outside of the invariant)
if length(t) > 1

    ind = find(list == 0);
    before = fun(t(end-1),x(:,end-1));
    after = fun(t(end),x(:,end));
    small = find(abs(before(ind)) < eps | abs(after(ind)) < eps);
    
    if any(before(ind) > 0 & after(ind) < 0) || ...
           (~isempty(small) && any(before(ind(small)) > after(ind(small))))
        nextloc = currentLoc; xJump = x(:,end); return;
    end
end

% loop over all active guards
for iActivatedGuard = 1:nActivatedGuards
    
    guard = list(index(iActivatedGuard));
    
    if guard ~= 0 && (~any(activatedGuards == guard))
        
        guardSet = loc.transition(guard).guard;
        
        % check whether only one halfspace has beed crossed or if 
        % the point is indeed inside the guard set
        if contains_(guardSet,x(:,end),'exact',1e-6,0,false,false)
            
            % next location and reset function
            nextloc = loc.transition(guard).target;
            xJump = evaluate(loc.transition(guard).reset,x(:,end));
            break;
            
        else
            
            % only one halfspace crossed -> continue simulation in
            % the current location
            nextloc = currentLoc;
            xJump = x(:,end);
        end
        
        activatedGuards = [activatedGuards,guard];
    end
end

% ------------------------------ END OF CODE ------------------------------
