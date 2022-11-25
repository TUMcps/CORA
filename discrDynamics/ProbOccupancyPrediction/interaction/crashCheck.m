function [crash] = crashCheck(carFstate,carLstate,carFinput,carLinput,T)
% crashCheck - checks if car F (follower) can crash into car L (leader) 
% if both cars accelerate with the given acceleration commands carFinput
% and carLinput for a time horizon T and an emergency brake with maxAcc 
% of both cars afterwards.
%
% Syntax:  
%    [crash] = crashCheck(carFstate,carLstate,carFinput,carLinput,T)
%
% Inputs:
%    carFstate - initial state of car F 
%    carLstate - initial state of car L
%    carFinput - input of car F
%    carLinput - input of car L
%    T - time interval for carFinput, carLinput
%
% Outputs:
%    crash - 1 if crash occurs and 0 otherwise
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      26-October-2007 
% Last update:  27-June-2008
%               12-October-2009
% Last revision:---

%------------- BEGIN CODE --------------

%combine initial states and inputs
x0=[carFstate(1); carLstate(1); carFstate(2); carLstate(2);];
u=[carFinput;carLinput];

%define event function from halfspace inequalities
eventOptions = odeset('Events',eventFcn());
%change other options
stepsizeOptions = odeset('MaxStep',0.2*T);
%generate overall options
options = odeset(eventOptions,stepsizeOptions);

%acceleration for $t\in[0,T]$
[t,x,te,xe,index] = ode45(getfcn(u),[0, T],x0,options);

if isempty(index)
    %set new initial states, inputs
    x0=x(end,:);
    u=[-1;-1];
    %simulate braking
    [t,x,te,xe,index] = ode45(getfcn(u),[0, 100*T],x0,options);
end

if index==1
    crash=1;
else
    crash=0;
end
    
end

%return function to be simulated
function [handle] = getfcn(u)
    function dxdt = f(t,x)
        dxdt = twoVehEid(t,x,u);
    end
    handle = @f;
end


%standard syntax for event functions as needed in Matlab simulations
function [handle] = eventFcn()
    function [value,isterminal,direction] = f(t,x)
        value(1)=x(2)-x(1);
        value(2)=x(3);
        % Always stop the integration when event detected
        isterminal = [1 1];   
        % direction of zero hitting; 0: both directions
        direction = [0 0];  
    end

    %return function handle
    handle = @f;

end

%------------- END OF CODE --------------