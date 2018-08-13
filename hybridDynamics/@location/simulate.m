function [t,x,loc,xJump] = simulate(obj,opt,tstart,tfinal,x0)
% simulate - simulates the system within a location, detects the guard set
% that is hit and computes the reset
%
% Syntax:
%    [t,x,loc,xJump] = simulate(obj,tstart,tfinal,x0)
%
% Inputs:
%    obj - location object
%    tstart - start time
%    tfinal - final time
%    x0 - initial state
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

%if there exists continuous dynamics
if ~isempty(obj.contDynamics)
    %define event function from halfspace inequalities
    eventOptions = odeset('Events',eventFcn(obj));
    %step-size options
    %stepsizeOptions = odeset('MaxStep',1e-2*(tstart-tfinal));
    %stepsizeOptions = odeset('MaxStep',0.2*(tstart-tfinal),'RelTol',1e-6,'AbsTol',1e-9);
    %generate overall options
    %options = odeset(eventOptions,stepsizeOptions);
    options = odeset(eventOptions);
    %simulate continuous dynamics
    [obj.contDynamics,t,x,index] = simulate(obj.contDynamics,opt,tstart,tfinal,x0,options);
    
    %time has not run out
    if ~isempty(index)
        %determine active guard
        [list] = indexList(obj);
        
        %         guard = list(index(end)); %take last event since other event were set when no termination was intended
        %         but events might be terminal and simultaneous....
        
        nActivatedGuards = length(index);
        activatedGuards = [];
        
        for iActivatedGuard = 1:nActivatedGuards
            
            guard = list(index(iActivatedGuard));
            
            if(guard == 0)
                %error('error in @location/simulate.m : invariant reached');
            elseif(~any(activatedGuards == guard))
                %check whether only one event function has been activated or if the
                %state is actually in the guard
                guardSet = get(obj.transition{guard},'guard');
                
                if in(guardSet, zonotope([x(end,:)', 1e3*eps*eye(length(x(end,:)))]))
                    %determine next location
                    loc=get(obj.transition{guard},'target');
                    %determine state reset
                    xJump=reset(obj.transition{guard},x(end,:));
                    break;
                else
                    loc = obj.id;
                    xJump = x(end,:);
                end
                activatedGuards = [activatedGuards,guard];
            end
        end
    else
        loc=[]; xJump=[];
    end
    
    %if there is no continuous dynamics
else
    loc=get(obj.transition{1},'target');
    t=tstart;
    x=x0';
    xJump=x0';
end

%------------- END OF CODE --------------
