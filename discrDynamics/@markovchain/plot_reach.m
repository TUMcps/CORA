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

% Author:       Matthias Althoff
% Written:      15-September-2006 
% Last update:  26-March-2008
%               24-July-2020
% Last revision: ---

%------------- BEGIN CODE --------------

%read objects
Obj = varargin{1};
HA = varargin{2};
R = varargin{3};
params = varargin{4};
actualSegmentNr = varargin{5};


%plot sample trajectories and reachable set 
h = figure;
set(gcf,'Units','normalized');
set(h,'position',[0.1,0.1,0.9,0.3]);
hold on

%start plotting
subplot(1,3,1); 
plot(Obj.field);
hold on
traj_plot(HA,R,params);

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
function traj_plot(HA,R,params)

% plot reachable set
plot(R,[1,2],'b','EdgeColor','b');

% plot initial set
plot(params.R0,[1,2],'w','Filled',true,'EdgeColor','k');

% settings for random simulation
simOpt.points = 30;        % number of initial points
simOpt.fracVert = 0.5;     % fraction of vertices initial set
simOpt.fracInpVert = 0.5;  % fraction of vertices input set
simOpt.inpChanges = 10;    % changes of input over time horizon  

% random simulation
simRes = simulateRandom(HA,params,simOpt); 

% plot simulated trajectories
plot(simRes,[1,2],'k');


%-------------------------------------------------------