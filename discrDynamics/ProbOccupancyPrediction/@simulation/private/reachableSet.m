function [indices,minPos,maxPos]=reachableSet(simOptions,markovChainSpec)
% reachableSet - computes the non-probabilistic reachable set of the
% vehicle
%
% Syntax:  
%    [indices,minPos,maxPos]=reachableSet(simOptions,markovChainSpec)
%
% Inputs:
%    simOptions - simulation options
%    markovChainSpec - Markov-Chain specifications
%
% Outputs:
%    indices - indices of reachable cells
%    minPos - minimum position of the vehicle
%    maxPos - maximum position of the vehicle
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       08-July-2008
% Last update:   28-August-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%get initial intervals
intervals=get(simOptions.initialStateSet,'intervals');

%obtain time vector
finalTime=simOptions.runs*markovChainSpec.timeStep;
steps=100;
deltaT=markovChainSpec.timeStep/steps;
t=0:deltaT:finalTime;

%simulate the braking trajectory
%conservative assumption: the vehicle has full negative acceleration
%potential, although it might start in a curve;
%the solution is obtained numerically, although it can be obtained
%analytically - however this approach is easier generalizable to less
%conservative computations

%set initial position
pos(1)=intervals(1,1);
vel(1)=intervals(2,1);
maxInput=-1;
acc=input2acceleration(maxInput,NaN,simOptions.type);

for i=1:length(t)
    %vehicle has not stopped yet
    if vel(i)>0
        pos(i+1)=pos(i)+vel(i)*deltaT;
        vel(i+1)=vel(i)+acc*deltaT;
    %vehicle has stopped    
    else
        pos(i+1)=pos(i);
        vel(i+1)=0;
    end
end

%get minimum positions, velocities
for i=1:simOptions.runs
    minPos(i)=pos((i-1)*steps+1);
    minVel(i)=vel((i-1)*steps+1);
end


%simulate the acceleration trajectory
%computation is based on the velocity profile
%use 'sliding mode approach'

%set initial position
pos(1)=intervals(1,2);
vel(1)=intervals(2,2);
maxInput=1;
minInput=-1;

for i=1:length(t)
    %get profile velocity
    profileVel=simOptions.profileHandle(pos(i),10);
    if vel(i)<profileVel
        acc=input2acceleration(maxInput,vel(i),simOptions.type);
    else
        acc=input2acceleration(minInput,vel(i),simOptions.type);;
    end
    pos(i+1)=pos(i)+vel(i)*deltaT;
    vel(i+1)=vel(i)+acc*deltaT;
end

%get maximum positions
for i=1:simOptions.runs
    maxPos(i)=pos((i-1)*steps+1);
    maxVel(i)=vel((i-1)*steps+1);
end

%obtain interval hulls of each time step and compute indices
for i=1:simOptions.runs
    IH=interval([minPos(i);minVel(i)],[maxPos(i);maxVel(i)]);
    ind=cellCandidates(simOptions.stateField,IH);
    [rows,cols,indices{i}]=find(ind); %pick only nonzero indices
end



%------------- END OF CODE --------------