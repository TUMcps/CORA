function plot_reach(varargin)
% plots - Generates 3 plots of a Markov chain:
% 1. Plot the continuous reachable set together with sample trajectories
% 2. Plot the reachable cells for the time point
% 3. Plot the reachable cells for the time interval
%
% Syntax:  
%    plot(Obj,HA,options,(actualSegmentNr))
%
% Inputs:
%    Obj - markovchain object
%    HA - hybrid automaton object
%    options - options struct
%    actualSegmentNr - number of the actual cell of the discretized state
%    space
%
% Outputs:
%    ---
%
% Example: 
%    ---
%
% Other m-files required: plotP, plot (for HA)
% Subfunctions: traj_plot
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 15-September-2006 
% Last update: 26-March-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%read objects
Obj=varargin{1};
HA=varargin{2};
options=varargin{3};
actualSegmentNr=varargin{4};


%plot sample trajectories and reachable set 
h=figure;
set(gcf,'Units','normalized');
set(h,'position',[0.1,0.1,0.9,0.3]);
hold on

%start plotting
subplot(1,3,1); 
plot(Obj.field);
hold on
traj_plot(HA,options);

%choose input that has been devoloped in the end
iInput=length(Obj.T.T);

T=Obj.T.T{iInput}; 
subplot(1,3,2);
plotP(Obj,T(:,actualSegmentNr+1),'k');
plot(Obj.field);
xlabel('x_1');
ylabel('x_2');

T=Obj.T.OT{iInput};
subplot(1,3,3); 
plotP(Obj,T(:,actualSegmentNr+1),'k');
plot(Obj.field);
xlabel('x_1');
ylabel('x_2');

%-------------------------------------------------------
%traj_plot: generates sample trajectories
function traj_plot(HA,options)

%plot reachable set
options.plotType = 'b';
plot(HA,'reachableSet',options);

%initial set
plot(options.R0,[1 2],'k-', 'lineWidth', 3);

%obtain random simulation results
for i=1:30
    %set initial state, input
    options.x0=randPoint(options.R0); %initial state for simulation
    options.uLoc{1}=randPoint(options.Uloc{1})+options.uLocTrans{1}; %input for simulation
    
    %simulate hybrid automaton
    HAsim{i}=simulate(HA,options);      
end

%plot simulation results      
for i=1:length(HAsim)
    plot(HAsim{i},'simulation',options);
end


%-------------------------------------------------------