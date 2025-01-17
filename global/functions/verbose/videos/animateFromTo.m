function handles = animateFromTo(vidObj,animationStruct,duration,varargin)
% animateFromTo - animates sets over the specified duration
%
% Syntax:
%    animateFromTo(vidObj,animationStruct,duration,frameRate)
%
% Inputs:
%    vidObj - VideoWriter
%    animationStruct - struct, specifying the animation; fields:
%       .fromS - contSet or numeric
%       .dims - numeric, dimensions to plot
%       .fromNVpairs - cell, name-value pairs for plotting
%       either:
%           .toS - contSet or numeric
%           .transformFun - function handle @(S,t) with time t \in [0,1]
%       .toNVpairs - (optional) cell, updated name-value pairs for plotting
%       .plotInBackground (optional) logical, whether to plot in background
%    duration - numeric positive scalar
%    varargin - name-value pairs
%       <'FromView',fromView> - numeric, [az,el] as required by view(az,el)
%       <'ToView',toView> - numeric, [az,el] as required by view(az,el)
%       <'FromEnlargeAxis',fromEnlargeAxis> - logical, whether to enlarge axis
%       <'ToEnlargeAxis',toEnlargeAxis> - logical, whether to enlarge axis
%       <'FreezeDurationBefore',freezeBefore> - numeric, duration to freeze animation beforehand
%       <'FreezeDurationAfter',freezeAfter> - numeric, duration to freeze animation afterward
%
% Outputs:
%    handles - cell, handles to the plotted graphics objects for t=1
%
% See also:
%    CORAvideo_snippets

% Authors:       Tobias Ladner
% Written:       10-December-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(3,Inf);
[fromView,toView,fromEnlargeAxis,toEnlargeAxis,freezeBefore,freezeAfter] = aux_readInput(vidObj,animationStruct,duration,varargin);
aux_parseInput(vidObj,animationStruct,duration,fromView,toView,fromEnlargeAxis,toEnlargeAxis)
animationStruct = aux_updateProperties(animationStruct);

% read out frame rate
FrameRate = vidObj.FrameRate;

% compute number of steps
nrOfSteps = duration * FrameRate;

% prepare animation
layer = nnSigmoidLayer();
progress = [0,layer.evaluate(linspace(-9,9,nrOfSteps-2)),1];

% get axis limits ---
[xLimFrom,yLimFrom,zLimFrom,xLimTo,yLimTo,zLimTo] = aux_gatherLimits(animationStruct,fromEnlargeAxis,toEnlargeAxis);

disp("- Start animating")

for i=1:nrOfSteps
    % ax = gca; fontSize = 24;
    % set(ax, "FontSize", fontSize)

    % get current progress of animation
    p = progress(i);

    fprintf('  Step %i/%i (p=%.4f) ..\n',i,nrOfSteps,p)
    
    % set limits
    xlim((1-p)*xLimFrom+p*xLimTo); xlabel('x_{(1)}')
    ylim((1-p)*yLimFrom+p*yLimTo); ylabel('x_{(2)}')
    zlim((1-p)*zLimFrom+p*zLimTo); zlabel('x_{(3)}')

    % init handle cell
    handles = cell(1,numel(animationStruct));

    for j=1:numel(animationStruct)
        % build in-between struct
        plotStruct = struct;
        plotStruct.S = animationStruct(j).transformFun(animationStruct(j).fromS,p);
        plotStruct.dims = animationStruct(j).dims;
        plotStruct.NVpairs = aux_mergeNVpairs(animationStruct(j),p);

        % plot set
        handles{j} = aux_plotSet(plotStruct);

        if animationStruct(j).plotInBackground
            uistack(handles{j},'bottom')
        end
    end

    % set view
    view_i = (1-p)*fromView+p*toView;
    view(view_i(1),view_i(2))

    % sort legend (to avoid weird jumps between animations)
    % childs = gca().Children;
    % [~,idx] = sort({childs.DisplayName});
    % legend(childs(idx))

    % finish current frame
    % legend('Location','northeastoutside')
    drawnow

    if i == 1
        % freeze video before animation
        writeFreezedFrames(vidObj,freezeBefore)
    end

    % write current frame
    writeVideo(vidObj, getframe(gcf));

    % delete all handles
    if i < nrOfSteps
        for j=1:numel(handles)
            delete(handles{j});
        end
    end
end

% freeze video after animation
writeFreezedFrames(vidObj,freezeAfter)

end


% Auxiliary functions -----------------------------------------------------

function [fromView,toView,fromEnlargeAxis,toEnlargeAxis,freezeBefore,freezeAfter] = aux_readInput(vidObj,animationStruct,duration,NVpairs)
    % read input

    % read name-value pairs
    checkNameValuePairs(NVpairs,{'FromView','ToView','FromEnlargeAxis','ToEnlargeAxis','FreezeDurationBefore','FreezeDurationAfter'});
    [NVpairs,fromView] = readNameValuePair(NVpairs,'FromView',{},[-35,30]);
    [NVpairs,toView] = readNameValuePair(NVpairs,'ToView','isnumeric',fromView);
    [NVpairs,fromEnlargeAxis] = readNameValuePair(NVpairs,'FromEnlargeAxis','islogical',false);
    [NVpairs,toEnlargeAxis] = readNameValuePair(NVpairs,'ToEnlargeAxis','islogical',false);
    [NVpairs,freezeBefore] = readNameValuePair(NVpairs,'FreezeDurationBefore','isnumeric',0);
    [NVpairs,freezeAfter] = readNameValuePair(NVpairs,'FreezeDurationAfter','isnumeric',0);
end

function aux_parseInput(vidObj,animationStruct,duration,fromView,toView,fromEnlargeAxis,toEnlargeAxis)
    % parse input

    % check classes
    inputArgsCheck({ ...
        {vidObj,'att','VideoWriter'}; ...
        {animationStruct,'att','struct'}; ...
        {duration,'att','numeric',{'positive','scalar'}}; ...
    })

    if ~isempty(animationStruct)
    
        % check struct fields
        requiredFields = {'fromS','dims','fromNVpairs'};
        for i=1:numel(requiredFields)
            field = requiredFields{i};
            if ~ismember(field,fieldnames(animationStruct))
                throw(CORAerror('CORA:wrongValue','second',sprintf('Struct ''animationStruct'' is missing the field ''%s''', field)))
            end
        end
        
        % check struct properties ---
        
        % .fromS
        func = @(animationStruct) isa(animationStruct.fromS,'contSet') || isnumeric(animationStruct.fromS) || isa(animationStruct.fromS,'reachSet');
        if ~all(arrayfun(func, animationStruct))
            throw(CORAerror('CORA:wrongValue','second','animationStruct.fromS must be of type contSet or numeric.'))
        end
    
        % .toS or transform given?
        func = @(animationStruct) (isfield(animationStruct,'toS') && ~isempty(animationStruct.toS)) || (isfield(animationStruct,'transformFun') && ~isempty(animationStruct.transformFun));
        if ~all(arrayfun(func, animationStruct))
            throw(CORAerror('CORA:wrongValue','second','animationStruct must have either fields .toS or .transformFun'))
        end
    
        % .toS
        func = @(animationStruct) ~isfield(animationStruct,'toS') || (isa(animationStruct.toS,'contSet') || isnumeric(animationStruct.toS) || isempty(animationStruct.toS));
        if ~all(arrayfun(func, animationStruct))
            throw(CORAerror('CORA:wrongValue','second','animationStruct.fromS must be of type contSet or numeric.'))
        end
    
        % .transform
        func = @(animationStruct) ~isfield(animationStruct,'transformFun') || isa(animationStruct.transformFun,'function_handle') || isempty(animationStruct.transformFun);
        if ~all(arrayfun(func, animationStruct))
            throw(CORAerror('CORA:wrongValue','second','animationStruct.transformFun must be of a function handle.'))
        end
        
        % .dims
        func = @(scalarStruct) isnumeric(scalarStruct.dims) && numel(scalarStruct.dims) <= 3;
        if ~all(arrayfun(func, animationStruct))
            throw(CORAerror('CORA:wrongValue','second','animationStruct.dims must specify at most three dimensions.'))
        end
        
        % .fromNVpairs
        func = @(animationStruct) iscell(animationStruct.fromNVpairs);
        if ~all(arrayfun(func, animationStruct))
            throw(CORAerror('CORA:wrongValue','second','animationStruct.fromNVpairs must be a cell array.'))
        end
    
        % .toNVpairs
        func = @(animationStruct) ~isfield(animationStruct,'toNVpairs') || isempty(animationStruct.toNVpairs) || iscell(animationStruct.toNVpairs);
        if ~all(arrayfun(func, animationStruct))
            throw(CORAerror('CORA:wrongValue','second','animationStruct.toNVpairs must be a cell array or empty.'))
        end
    
        % .plotInBackground
        func = @(animationStruct) ~isfield(animationStruct,'plotInBackground') || isempty(animationStruct.plotInBackground) || islogical(animationStruct.plotInBackground);
        if ~all(arrayfun(func, animationStruct))
            throw(CORAerror('CORA:wrongValue','second','animationStruct.plotInBackground must be a logical or empty.'))
        end
    end

    % check other properties ---
    
    % view
    func = @(view) isnumeric(view) && numel(view) == 2;
    if ~func(fromView) || ~func(toView)
        throw(CORAerror('CORA:wrongValue','name-value pair fromView/toView','View must have two numeric entries [az,el].'))
    end

    % enlargeAxis
    func = @(enlargeAxis) islogical(enlargeAxis) && isscalar(enlargeAxis);
    if ~func(fromEnlargeAxis) || ~func(toEnlargeAxis)
        throw(CORAerror('CORA:wrongValue','name-value pair fromEnlargeAxis/toEnlargeAxis','DoEnlargeAxis must be a logical scalar.'))
    end

end

function animationStruct = aux_updateProperties(animationStruct)

    % transform every .toS into .transformFun, make sure .toS is present
    for i=1:numel(animationStruct)
        if isfield(animationStruct(i), 'toS') && ~isempty(animationStruct(i).toS)
            % convex combination in between .fromS and .toS
            animationStruct(i).transformFun = @(S,t) (1-t)*S + t*animationStruct(i).toS;
        else
            % compute final set .toS
            animationStruct(i).toS = animationStruct(i).transformFun(animationStruct(i).fromS,1);
        end
    end

    % check if .toNVpairs are present
    for i=1:numel(animationStruct)
        if ~isfield(animationStruct(i),'toNVpairs') || isempty(animationStruct(i).toNVpairs)
            animationStruct(i).toNVpairs = animationStruct(i).fromNVpairs;
        end
    end

    % check if .plotInBackground is present
    for i=1:numel(animationStruct)
        if ~isfield(animationStruct(i),'plotInBackground') || isempty(animationStruct(i).plotInBackground)
            animationStruct(i).plotInBackground = false;
        end
    end

end

function [xLimFrom,yLimFrom,zLimFrom,xLimTo,yLimTo,zLimTo] = aux_gatherLimits(animationStruct,fromEnlargeAxis,toEnlargeAxis)
    
    if strcmp(xlim('mode'),'manual')
        % read from axis limits directly to avoid jumps
        xLimFrom=xlim; yLimFrom=ylim; zLimFrom=zlim;
        
        % continue with 'to'
        wStart = 2;    
    else
        % also determine 'from' axis
        wStart = 1;
    end

    % init cell for handles
    handles = cell(2,numel(animationStruct));

    whichS = {'fromS','toS'};
    for w = wStart:2

        if w == 2
            % set mode to auto for final transition
            xlim('auto'); ylim('auto'); zlim('auto');
        end

        % plot all set
        for i=1:numel(animationStruct)
            % build default plot struct
            plotStruct = struct;
            plotStruct.S = animationStruct(i).(whichS{w});
            plotStruct.dims = animationStruct(i).dims;
            plotStruct.NVpairs = animationStruct(i).fromNVpairs;
    
            % plot
            handles{w,i} = aux_plotSet(plotStruct);
        end
    
        % read limits
        if w == 1
            if fromEnlargeAxis
                % enlarge axis
                enlargeAxis();
            end
            % read axis limits
            xLimFrom=xlim; yLimFrom=ylim; zLimFrom=zlim;
        else
            if toEnlargeAxis
                % enlarge axis
                enlargeAxis();
            end
            % read axis limits
            xLimTo=xlim; yLimTo=ylim; zLimTo=zlim;
        end
   
        % delete handles
        for i=1:numel(handles)
            delete(handles{i})
        end
    end
end

function han = aux_plotSet(plotStruct)
    % plots a single set

    % read parameters
    S = plotStruct.S;
    dims = plotStruct.dims;
    NVpairs = plotStruct.NVpairs;
    
    % plot
    if isnumeric(S)
        han = plotPoints(S,dims,NVpairs{:});
    else
        han = plot(S,dims,NVpairs{:});
    end
end

function NVpairs = aux_mergeNVpairs(animationStruct,delta)

    % read NVpairs
    NVpairs = animationStruct.fromNVpairs; % default name-value pairs
    toNVpairs = animationStruct.toNVpairs; % name-value pairs to update

    if isempty(toNVpairs)
        return
    end

    % following arguments need merging
    requireMerging = {'Color','EdgeColor','FaceColor','FaceAlpha','EdgeAlpha'};
    for i=1:numel(requireMerging)
        name = requireMerging{i};

        % read value
        [NVpairs,valueFrom] = readNameValuePair(NVpairs,name);
        [toNVpairs,valueTo] = readNameValuePair(toNVpairs,name);

        % update 'from' if necessary
        if ~isempty(valueTo) && ~isempty(valueFrom)
            NVpairs = [NVpairs {name,(1-delta)*valueFrom + delta*valueTo}];
        elseif ~isempty(valueFrom)
            NVpairs = [NVpairs {name,valueFrom}];
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
