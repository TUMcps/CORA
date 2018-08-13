function [obj] = simulate_tracking(obj,options)
% simulate - simulates a hybrid automaton for a tracking task
%
% Syntax:  
%    [obj] = simulate_tracking(obj,options)
%
% Inputs:
%    obj - hybrid automaton object
%    options - simulation options
%
% Outputs:
%    obj - hybrid automaton object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      21-September-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%load data from options
tFinal=options.tFinal; %final time
finalLoc=options.finalLoc; %final location

%initialize variables, intermediate state, 
loc=options.startLoc; %actual location
xInter=options.x0; %intermediate state at transitions
t=[]; %time vector
x=[]; %state vector


%prepare successive simulations
inputChanges=ceil(tFinal/options.timeStepLoc{loc});
tFinal = options.tStart;

for iChange = 1:inputChanges
    %reset options
    tInter=tFinal;
    tFinal=tFinal+options.timeStepLoc{loc};
    options.u = options.uLocTrans{loc}(:,iChange);

    %while final time not reached
    while (tInter<tFinal) && (~isempty(loc)) && (loc~=finalLoc) 
        %simulate within the actual location
        [tNew,xNew,newLoc,xInter]=simulate(obj.location{loc},options,tInter,tFinal,xInter);
        %new intermediate time is last simulation time
        tInter=tNew(end);
        
        if ~isempty(newLoc)
            loc = newLoc;
            disp('new location');
        end

        %store results
        t((end+1):(end+length(tNew))) = tNew;
        x((end+1):(end+length(tNew(:,1))),:) = xNew;
        xInter = xNew(end,:);
    end
end

%save results to object structure
obj.result.simulation.t=t;
obj.result.simulation.x=x;
obj.result.simulation.location=[];


%------------- END OF CODE --------------