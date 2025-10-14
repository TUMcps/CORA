function res = recordReachableSet(vidObj,fig,varargin)
% recordReachableSet - records a video of the computed reachable set
%
% Syntax:
%    res = recordReachableSet(vidObj,varargin)
%
% Inputs:
%    vidObj - VideoWriter
%    fig - figure
%    varargin - name-value-pairs
%       'ReachSet' - reachSet object
%       'Trajectory' - trajectory object
%       'RefTrajectory' - trajectory object
%       'Dimensions' - dimensions to plot (1 for plotOverTime, 2 for plot)
%       'RefDimensions' - dimensions to plot reference trajectory
%       'Specification' - specification object
%       'UnifyTotalSets' - total number of sets for unify
%       'TotalDuration' - total duration of video
%       'FreezeDuration' - duration of freezed animation at the end
%       'ReachSets' - cell array of {<reachSet>,<display-name>,<color>}
%       'PlotMethodTrajectory' - one of 'Time', 'Percent', 'AllAtOnce'
%
% Outputs:
%    res - logical
%
% See also:
%    CORAvideo_snippets

% Authors:       Tobias Ladner
% Written:       05-April-2024
% Last update:   19-April-2024 (TL, added reference trajectory)
%                11-December-2024 (TL, split setup and recording)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp('Recording reachable set:')

% parse input ---

[Rs,traj,refTrajectory,dims,refDims,spec, ...
    Unify,UnifyTotalSets,MaxTimeRs,TotalDuration,FreezeDuration,PlotMethodTrajectory] = aux_parseInput(varargin);
simName = 'Simulations';
refName = 'Ref. trajectory';

% setup ---

disp("- Setup")

% labels
if isscalar(dims)
    XLabel = 'Time';
else
    XLabel = sprintf('x_{(%i)}',dims(1));
end
xlabel(XLabel);
ylabel(sprintf('x_{(%i)}',dims(end)))

% enlarge axis in case not set already, ensures reachable set is doesn't
% change axis limits, which does not look nice on the video
if ~isempty([traj.x])
    aux_plot(traj,dims,'DisplayName',simName,'Color',CORAcolor('CORA:simulations'));
    aux_enlargeAxis(1.2,dims);
    if ~strcmp(PlotMethodTrajectory,'all')
        aux_deletegraphics(simName);
    end
end

% plot initial set and simulations for axis scaling
specName = 'Goal set';
aux_plot(spec,dims,'DisplayName',specName);
if ~isscalar(dims)
    % do not plot for plotOverTime
    aux_plot(Rs{1,1}.R0,dims,'DisplayName','Initial set','Color',CORAcolor("CORA:initialSet"));
end

% record video ---

disp("- Start recording")

% first frame
writeVideo(vidObj, getframe(gcf));

% get time per step
FrameRate = vidObj.FrameRate;
nrFrames = TotalDuration * FrameRate;
nrFramesFreeze = FreezeDuration * FrameRate;
nrFramesAnimation = nrFrames - nrFramesFreeze;
deltat = (MaxTimeRs/(TotalDuration-FreezeDuration)) / FrameRate;

% init default reachSet name-value pairs
reachSetNVpairs = {'Unify',Unify,'UnifyTotalSets',UnifyTotalSets, ...
            'PlotBackground',true};
tStart = 0;
tFinalUpToFrame = 0;
Rs_to_i = cell(1,size(Rs,1));

% plot reachable set
for frame=1:nrFramesAnimation

    % init
    t_i = deltat * (frame-1);
    t_i1 = deltat * frame;
    tau_0i1 = interval(tStart,t_i1);
    fprintf('  [%.2f, %.2f] / %.2f ..\n', t_i, t_i1, MaxTimeRs);
    title(sprintf('$t=%.2fs$',t_i1),'Interpreter','latex')

    % plot reachable set
    timerVal = tic;
    for r=1:size(Rs,1)
        R_all = Rs{r,1};

        % get all reachable sets of this frame
        R_to_i = find(R_all,'time',tau_0i1);

        % reach properties
        R_name = Rs{r,2};
        R_color = Rs{r,3};

        % remove old reachable sets
        aux_deletegraphics(R_name);

        % plot reachable sets
        aux_plot(R_to_i,dims,'DisplayName',R_name,'Color',R_color,reachSetNVpairs{:});
        aux_moveToBackground(R_name)

        % update max time
        tFinalUpToFrame = max(tFinalUpToFrame,query(R_to_i,'tFinal'));

        % save for next round
        Rs_to_i{r} = R_to_i;
    end
    plotTime = toc(timerVal);

    % if time to plot reachable sets took too long, 
    % keep current reachable sets and only delete subsequent
    if plotTime > size(Rs,1) * 0.1 % s
        for r=1:size(Rs,1)
            % add underscore to distinguish new sets by name
            Rs{r,2} = [Rs{r,2} '_'];
        end
        reachSetNVpairs = [reachSetNVpairs, {'HandleVisibility','off'}];
        tStart = min(t_i, tFinalUpToFrame);
    end

    % plot simulations (delete all previous simulations)
    if ~isempty([traj.x])
        if strcmp(PlotMethodTrajectory, 'time')
            % plot simulations at the same time as reachable set
            aux_deletegraphics(simName);
            tau_0i1 = interval(0,min(t_i1,tFinalUpToFrame));
            aux_plot(find(traj,'time',tau_0i1),dims, ...
                'DisplayName',simName,'Color',CORAcolor('CORA:simulations'));
    
        elseif strcmp(PlotMethodTrajectory, 'percent')
            % iteratively plot new simulations based on the current progress
            % (similar to website)
            aux_deletegraphics(simName);
            percent = ceil(frame/nrFrames * numel(traj));
            aux_plot(traj(1:percent),dims, ...
                'DisplayName',simName,'Color',CORAcolor('CORA:simulations'));
        end
    end
    
    % plot reference trajectory
    if ~isempty(refTrajectory)
        if strcmp(PlotMethodTrajectory, 'time')
            % plot reference trajectory at the same time as reachable set
            aux_deletegraphics(refName);
            tau_0i1 = interval(0,min(t_i1,tFinalUpToFrame));
            aux_plot(find(refTrajectory,'time',tau_0i1),refDims, ...
                'DisplayName',refName,'Color',[1 0 0],'LineWidth',2);
    
        elseif strcmp(PlotMethodTrajectory, 'percent')
            % always plot single given reference trajactory
            % (similar to website)
            aux_deletegraphics(refName);
            aux_plot(refTrajectory,refDims, ...
                'DisplayName',refName,'Color',[1 0 0],'LineWidth',2);
        end
    end

    % move specification back
    aux_moveToBackground(specName)

    % write frame
    writeVideo(vidObj, getframe(gcf));
end

% write freezed frames
disp('- Freezed frames')
writeFreezedFrames(vidObj,FreezeDuration);

% completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function [Rs,traj,refTrajectory,dims,refDims,spec,Unify,UnifyTotalSets,MaxTimeRs,TotalDuration,FreezeDuration,PlotMethodTrajectory] = aux_parseInput(NVpairs)
    
    % read input args
    [NVpairs,R] = readNameValuePair(NVpairs,'ReachSet');
    [NVpairs,Rs] = readNameValuePair(NVpairs,'ReachSets');
    [NVpairs,traj] = readNameValuePair(NVpairs,'Trajectory');
    [NVpairs,refTrajectory] = readNameValuePair(NVpairs,'RefTrajectory');
    [NVpairs,dims] = readNameValuePair(NVpairs,'Dimensions');
    [NVpairs,refDims] = readNameValuePair(NVpairs,'RefDimensions');
    [NVpairs,spec] = readNameValuePair(NVpairs,'Specification');
    [NVpairs,Unify] = readNameValuePair(NVpairs,'Unify');
    [NVpairs,UnifyTotalSets] = readNameValuePair(NVpairs,'UnifyTotalSets');
    [NVpairs,FreezeDuration] = readNameValuePair(NVpairs,'FreezeDuration');
    [NVpairs,TotalDuration] = readNameValuePair(NVpairs,'TotalDuration');
    [NVpairs,PlotMethodTrajectory] = readNameValuePair(NVpairs,'PlotMethodTrajectory');

    % check must haves
    if isempty(R) && isempty(Rs)
        throw(CORAerror('CORA:wrongValue','name-value pair ''ReachSet''','missing mandatory value'))
    end
    if isempty(traj)
        throw(CORAerror('CORA:wrongValue','name-value pair ''Trajectory''','missing mandatory value'))
    end
    if isempty(dims)
        throw(CORAerror('CORA:wrongValue','name-value pair ''Dimensions''','missing mandatory value'))
    end

    % add R to Rs
    if ~isempty(R)
        Rs = [{R,'Reachable set',CORAcolor("CORA:reachSet")};Rs];
    end

    % flip Rs for intuitive understanding due to background plotting
    Rs = flipud(Rs);

    % check reference trajectory
    if ~isempty(refTrajectory) && ~isa(refTrajectory,'trajectory')
        throw(CORAerror('CORA:wrongValue','name-value pair ''RefTrajectory''','trajectory'))
    end

    % set missing values ---

    if isempty(refDims)
        refDims = dims;
    end

    % unify
    if isempty(Unify)
        Unify = true;
    end
    if isempty(UnifyTotalSets)
        UnifyTotalSets = 1;
    end

    % maximum time of any reachable set
    MaxTimeRs = 0;
    for i=1:size(Rs,1)
        MaxTimeRs = max(MaxTimeRs,query(Rs{i,1},'tFinal'));
    end

    % video duration
    if isempty(FreezeDuration)
        FreezeDuration = 1;
    end
    if isempty(TotalDuration)
        TotalDuration = MaxTimeRs + FreezeDuration;
    end

    % ploting method for simulations
    if isempty(PlotMethodTrajectory)
        PlotMethodTrajectory = 'time';
    end
end

function han = aux_plot(obj, dims, varargin)
    % only plot if not empty
    if ~isempty(obj)
        % plotOverTime vs. standard plot 
        if isscalar(dims)
            han = plotOverTime(obj,dims,varargin{:});
        else
            han = plot(obj,dims,varargin{:});
        end
    end
end

function aux_enlargeAxis(factor,dims)
    
    % get axis limits
    xLim = xlim();
    yLim = ylim();
    
    % enlarge viewbox
    I = interval([xLim(1);yLim(1)],[xLim(2);yLim(2)]);
    I = enlarge(I,factor);
    
    % set axis
    if strcmp(xlim('mode'),'auto')
        if isscalar(dims)
            % keep time axis
            xlim(xLim)
        else
            xlim([I.inf(1) I.sup(1)])
        end
    end
    if strcmp(ylim('mode'),'auto')
        ylim([I.inf(2) I.sup(2)])
    end
end

function aux_deletegraphics(name)
    children = allchild(gca);
    idx = strcmp(get(children,'DisplayName'),name);
    delete(children(idx));
end

function aux_moveToBackground(name)
    % remove added underscores
    while endsWith(name,'_')
        name = name(1:end-1);
    end

    % move all children with that name back
    children = allchild(gca);
    idx = startsWith(get(children,'DisplayName'),name);
    uistack(children(idx),'bottom');
end

% ------------------------------ END OF CODE ------------------------------
