function res = recordCORAvideo(filename,varargin)
% recordCORAvideo - records a video of the computed reachable set
%
% Syntax:
%    res = recordCORAvideo(filename,varargin)
%
% Inputs:
%    filename - char, video file name (with extension *.mp4)
%    varargin - name-value-pairs
%       'ReachSet' - reachSet object
%       'SimResult' - simResult object
%       'RefTrajectory' - simResult object
%       'Dimensions' - dimensions to plot (1 for plotOverTime, 2 for plot)
%       'RefDimensions' - dimensions to plot reference trajectory
%       'Specification' - specification object
%       'XLabel' - xlabel
%       'XLim' - xlim
%       'YLabel' - ylabel
%       'YLim' - ylim
%       'Title' - title of video
%       'Description' - description (also controls layout)
%       'LegendLocation' - location of legend
%       'UnifyTotalSets' - total number of sets for unify
%       'FrameRate' - frame rate of video
%       'TotalDuration' - total duration of video
%       'FreezeDuration' - duration of freezed animation at the end
%       'ReachSets' - cell array of {<reachSet>,<display-name>,<color>}
%       'PlotMethodSimResult' - one of 'Time', 'Percent', 'AllAtOnce'
%
% Outputs:
%    res - logical
%
% Example
%    % gather R and simRes
%    % ...
%    recordCORAvideo('example.mp4','ReachSet',R,'SimResult',simRes,'Dimensions',[1]);

% Authors:       Tobias Ladner
% Written:       05-April-2024
% Last update:   19-April-2024 (TL, added reference trajectory)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp('Recording CORA video:')

% parse input ---

[Rs,simRes,refTrajectory,dims,refDims,spec, ...
    XLabel,XLim,YLabel,YLim,Title,Description,LegendLocation, ...
    Unify,UnifyTotalSets, ...
    FrameRate,MaxTimeRs,TotalDuration,FreezeDuration,PlotMethodSimResult] = aux_parseInput(varargin);
simName = 'Simulations';
refName = 'Ref. trajectory';

% setup ---

disp("- Setup")

% set up figure
fig = aux_setUpFigure(XLabel,XLim,YLabel,YLim,Title,Description,LegendLocation);

% set up video
vidObj = VideoWriter(filename, 'MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = FrameRate;
open(vidObj);

% enlarge axis in case not set already, ensures reachable set is doesn't
% change axis limits, which does not look nice on the video
aux_plot(simRes,dims,'DisplayName',simName,'Color',CORAcolor('CORA:simulations'));
aux_enlargeAxis(1.2,dims);
if ~strcmp(PlotMethodSimResult,'all')
    aux_deletegraphics(simName);
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
    fprintf('  [%.2f, %.2f] / %.2f \n', t_i, t_i1, MaxTimeRs);

    % plot reachable set
    tic
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
    plotTime = toc;

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
    if strcmp(PlotMethodSimResult, 'time')
        % plot simulations at the same time as reachable set
        aux_deletegraphics(simName);
        tau_0i1 = interval(0,min(t_i1,tFinalUpToFrame));
        aux_plot(find(simRes,'time',tau_0i1),dims, ...
            'DisplayName',simName,'Color',CORAcolor('CORA:simulations'));

    elseif strcmp(PlotMethodSimResult, 'percent')
        % iteratively plot new simulations based on the current progress
        % (similar to website)
        aux_deletegraphics(simName);
        percent = ceil(frame/nrFrames * numel(simRes));
        aux_plot(simRes(1:percent),dims, ...
            'DisplayName',simName,'Color',CORAcolor('CORA:simulations'));
    end

    % plot reference trajectory
    if ~isempty(refTrajectory)
        if strcmp(PlotMethodSimResult, 'time')
            % plot reference trajectory at the same time as reachable set
            aux_deletegraphics(refName);
            tau_0i1 = interval(0,min(t_i1,tFinalUpToFrame));
            aux_plot(find(refTrajectory,'time',tau_0i1),refDims, ...
                'DisplayName',refName,'Color',[1 0 0],'LineWidth',2);
    
        elseif strcmp(PlotMethodSimResult, 'percent')
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
writeVideo(vidObj, repmat(getframe(gcf),nrFramesFreeze,1));

% close
close(fig)
close(vidObj);

% res
fprintf('Done recording. Video saved as ''%s''.\n\n', filename)
res = true;

end


% Auxiliary functions -----------------------------------------------------

function [Rs,simRes,refTrajectory,dims,refDims,spec,XLabel,XLim,YLabel,YLim,Title,Description,LegendLocation,Unify,UnifyTotalSets,FrameRate,MaxTimeRs,TotalDuration,FreezeDuration,PlotMethodSimResult] = aux_parseInput(NVpairs)
    
    % read input args
    [NVpairs,R] = readNameValuePair(NVpairs,'ReachSet');
    [NVpairs,Rs] = readNameValuePair(NVpairs,'ReachSets');
    [NVpairs,simRes] = readNameValuePair(NVpairs,'SimResult');
    [NVpairs,refTrajectory] = readNameValuePair(NVpairs,'RefTrajectory');
    [NVpairs,dims] = readNameValuePair(NVpairs,'Dimensions');
    [NVpairs,refDims] = readNameValuePair(NVpairs,'RefDimensions');
    [NVpairs,spec] = readNameValuePair(NVpairs,'Specification');
    [NVpairs,XLabel] = readNameValuePair(NVpairs,'XLabel');
    [NVpairs,XLim] = readNameValuePair(NVpairs,'XLim');
    [NVpairs,YLabel] = readNameValuePair(NVpairs,'YLabel');
    [NVpairs,YLim] = readNameValuePair(NVpairs,'YLim');
    [NVpairs,Title] = readNameValuePair(NVpairs,'Title');
    [NVpairs,Description] = readNameValuePair(NVpairs,'Description');
    [NVpairs,LegendLocation] = readNameValuePair(NVpairs,'LegendLocation');
    [NVpairs,Unify] = readNameValuePair(NVpairs,'Unify');
    [NVpairs,UnifyTotalSets] = readNameValuePair(NVpairs,'UnifyTotalSets');
    [NVpairs,FrameRate] = readNameValuePair(NVpairs,'FrameRate');
    [NVpairs,FreezeDuration] = readNameValuePair(NVpairs,'FreezeDuration');
    [NVpairs,TotalDuration] = readNameValuePair(NVpairs,'TotalDuration');
    [NVpairs,PlotMethodSimResult] = readNameValuePair(NVpairs,'PlotMethodSimResult');

    % check must haves
    if isempty(R) && isempty(Rs)
        throw(CORAerror('CORA:wrongValue','name-value pair ''ReachSet''','missing mandatory value'))
    end
    if isempty(simRes)
        throw(CORAerror('CORA:wrongValue','name-value pair ''SimResult''','missing mandatory value'))
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
    if ~isempty(refTrajectory) && ~isa(refTrajectory,'simResult')
        throw(CORAerror('CORA:wrongValue','name-value pair ''RefTrajectory''','simResult'))
    end

    % set missing values ---

    if isempty(refDims)
        refDims = dims;
    end

    % labels
    if isempty(XLabel)
        if isscalar(dims)
            XLabel = 'Time';
        else
            XLabel = sprintf('x_{(%i)}',dims(1));
        end
    end
    if isempty(YLabel)
        YLabel = sprintf('x_{(%i)}',dims(end));
    end

    % legend
    if isempty(LegendLocation)
        LegendLocation = 'northeast';
    end

    % unify
    if isempty(Unify)
        Unify = true;
    end
    if isempty(UnifyTotalSets)
        UnifyTotalSets = 1;
    end

    % frame rate
    if isempty(FrameRate)
        FrameRate = 30;
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
    if isempty(PlotMethodSimResult)
        PlotMethodSimResult = 'time';
    end
end

function fig = aux_setUpFigure(XLabel,XLim,YLabel,YLim,Title,Description,LegendLocation)
    
    % set up figure
    fig = figure; 
    set(fig,'WindowStyle','normal');
    width = 2160;
    height = 2160;

    % background color white
    backgroundColor = [1 1 1];
    set(gcf, "Color", backgroundColor);

    % font size
    fontSize = 36;
    fontSizeTitle = ceil(fontSize * 1.33);
    fontSizeSmall = ceil(fontSize / 1.33);

    if isempty(Description)
        % just show plot
        aux_setFigureSize(fig, width, height)

        % and write title above plot
        hold on; box on;
        titleHandle = title(Title);
        titleHandle.FontSize = fontSizeTitle;
    else
        % make space for title and description
        width = 3840;
        aux_setFigureSize(fig, width, height)

        % plot title
        dim = [.025 .8 .45 .1]; 
        annotation('textbox',dim,'String',Title,'EdgeColor','none','FontWeight','bold','FontSize',fontSizeTitle)
        
        % plot description
        dim = [.025 .3 .45 .5]; 
        annotation('textbox',dim,'String',compose(Description),'EdgeColor','none','FontSize',fontSize)

        % show image
        img = imread('coraLogo.png','BackgroundColor',[1 1 1]);
        s = size(img);
        imgwidth = 0.075;
        axCORA = axes( ...
            'Units','pixels', ...
            'OuterPosition',[width*0.025 height*0.05 width*imgwidth s(1)/s(2)*imgwidth*width+fontSizeSmall/0.75] ...
        );
        image(axCORA,img)
        axCORA.Toolbar.Visible = 'off';
        disableDefaultInteractivity(axCORA)
        set(axCORA,'xtick',[])
        set(axCORA,'ytick',[])
        set(axCORA,'XColor',backgroundColor,'YColor',backgroundColor,'TickDir','out')
        website = xlabel(axCORA, 'https://cora.in.tum.de');
        website.Color = 'k';
        website.FontSize = fontSizeSmall;

        % create subplot on the right
        subplot(1,2,2); hold on; box on;
    end

    fprintf('- Quality: %ix%i\n', width, height)

    % set up axis to plot reachable set
    ax = gca();
    disableDefaultInteractivity(ax)
    ax.Toolbar.Visible = 'off';
    set(ax, "FontSize", fontSize)
    
    xlabel(XLabel); 
    if ~isempty(XLim)
        xlim(XLim)
    end
    ylabel(YLabel); 
    if ~isempty(YLim)
        ylim(YLim)
    end
    legend('Location',LegendLocation)
end

function aux_setFigureSize(fig, width, height)
    screenSize = get(0, 'ScreenSize');
    left = (screenSize(3) - width) / 2;
    bottom = 0;
    set(fig, 'OuterPosition', [left, bottom, width, height]);
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
