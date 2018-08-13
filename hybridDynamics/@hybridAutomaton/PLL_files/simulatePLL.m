function [obj] = simulatePLL(obj,options)
% simulate - simulates a hybrid automaton
%
% Syntax:  
%    [obj] = simulate(obj,tstart,tfinal,startlocation,y0)
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
% Written:      20-January-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%load data from options
tFinal=options.tFinal; %final time

%initialize variables, intermediate state, 
tInter=options.tStart; %intermediate time at transitions
loc=options.startLoc; %actual location
xInter=options.x0; %intermediate state at transitions
t=[]; %time vector
x=[]; %state vector
x_cycle=xInter'; %state vector after one cycle
location=[];
count=1; %transition counter

%projection matrix
P = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 -1 0];

%generate random inputs for each location
options.uLoc{1} = options.Uloc{1};
options.uLoc{2} = randPoint(options.Uloc{2});
options.uLoc{3} = randPoint(options.Uloc{3});
uTmp = options.Uloc{1};
uTmp(1:3) = options.uLoc{2}(1:3) + options.uLoc{3}(1:3);
options.uLoc{4} = uTmp;

%convert state into dummy interval
IH_state = interval(P*xInter,P*xInter);

%cycle number when locking region is reached
lockCycle=[];

%while final time not reached
while length(x_cycle(:,1)) < options.simCycles
    %store locations
    location(count)=loc;
    
    options.u = options.uLoc{loc};

    %simulate within the actual location
    [tNew,xNew,loc,xInter]=simulate(obj.location{loc},options,tInter,tFinal,xInter);
    %new intermediate time is last simulation time
    tInter=tNew(end);
    
    %store results
    t{count}=tNew;
    x{count}=xNew;
    
    %increase counter for transitions
    count=count+1;
    
    %store results at the end of location 1
    if loc==2 && (xInter(4)<=xInter(5))
        x_cycle(end+1,:) = xInter';
    end
    if loc==4 && (xInter(4)>=xInter(5))
        x_cycle(end+1,:) = xInter';
    end
    %convert state into dummy interval
    IH_state = interval(P*xInter,P*xInter);
    
    %check if locking condition is fulfilled
    if isempty(lockCycle)
        if IH_state<=options.IHstart
            lockCycle = length(x_cycle(:,1));
        end
    end
end


%save results to object structure
obj.result.simulation.t=t;
obj.result.simulation.x=x;
obj.result.simulation.x_cycle=x_cycle;
obj.result.simulation.lockCycle=lockCycle;
obj.result.simulation.location=location;


%------------- END OF CODE --------------