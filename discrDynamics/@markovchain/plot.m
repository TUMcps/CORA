function plot(varargin)
% plot - Generates 3 plots of a Markov chain:
% 1. Plot the sample trajectories
% 2. Plot the reachable cells for the time point
% 3. Plot the reachable cells for the time interval
%
% Syntax:  
%    plot(Obj,HA,options,(actualSegmentNr))
%
% Inputs:
%    Obj - markovchain object
%    trajectories - sample trajectories
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
%               16-June-2009
% Last revision: ---

%------------- BEGIN CODE --------------

%read objects
Obj=varargin{1};
trajectories=varargin{2};
options=varargin{3};
actualSegmentNr=varargin{4};

%plot sample trajectories and probability distribution
h=figure;
set(gcf,'Units','normalized');
set(h,'position',[0.1,0.1,0.9,0.3]);
hold on

% plot sample trajectories
subplot(1,3,1);
plot(Obj.field);
traj_plot(trajectories,options);

%choose input that has been devoloped in the end
iInput=length(Obj.T.T);

% plot distribution for points in time
T=Obj.T.T{iInput}; 
subplot(1,3,2);
plotP(Obj,T(:,actualSegmentNr+1),'k');
plot(Obj.field);
xlabel('x_1');
ylabel('x_2');

% plot distribution for time intervals
T=Obj.T.OT{iInput};
subplot(1,3,3);
plotP(Obj,T(:,actualSegmentNr+1),'k');
plot(Obj.field);
xlabel('x_1');
ylabel('x_2');

%-------------------------------------------------------
%traj_plot: generates sample trajectories
function traj_plot(trajectories,options)

%plot initial set
plot(options.R0,[1 2],'w','EdgeColor','k','Filled',true); 

hold on

%plot trajectories
for initIndx=1:length(trajectories(:,1))
    for inputIndx=1:length(trajectories(1,:))
        for iLoc=1:length(trajectories{initIndx,inputIndx})
            val=trajectories{initIndx,inputIndx}{iLoc};
            plot(val(:,1),val(:,2),'Color',0.6*[1 1 1]);
        end
    end
end

xlabel('x_1');
ylabel('x_2');
%-------------------------------------------------------