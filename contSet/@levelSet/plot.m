function han = plot(ls,varargin)
% plot - plots a projection of a level set
%
% Syntax:
%    han = plot(ls)
%    han = plot(ls,dims)
%    han = plot(ls,dims,type)
%
% Inputs:
%    ls - levelSet object
%    dims - (optional) dimensions for projection
%           (assumption: other entries of the normal vector are zeros)
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%               - 'Splits': number of splits to plot sets with '<=','<'
%               - 'PlotMethod': one of {'outer','inner'}, default: 'outer'
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    syms x y
%    eq = x^2 + y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%    
%    figure; hold on; box on; xlim([-3,3]); ylim([-3,3]);
%    plot(ls,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       19-July-2019
% Last update:   22-May-2023 (TL, speed up plotting of '<=' levelSets)
%                26-July-2023 (TL, getUnboundedAxisLimits)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input arguments
[ls,dims,NVpairs,splits,plotMethod] = aux_parseInput(ls,varargin{:});

% 2. preprocess
[ls,dims] = aux_preprocess(ls,dims);

% 3. plot n-dimensional set
han = aux_plotNd(ls,dims,NVpairs,splits,plotMethod);

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [ls,dims,NVpairs,splits,plotMethod] = aux_parseInput(ls,varargin)
    % parse input arguments
    dims = setDefaultValues({[1,2]},varargin);

    % check input arguments
    inputArgsCheck({{ls,'att','levelSet'};
                    {dims,'att','numeric',{'nonnan','integer','vector','positive'}}});

    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end

    % read out plot options and additional name-value pairs
    NVpairs = readPlotOptions(varargin(2:end));
    % 'Splits' given?
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits');
    if isempty(splits)
        if length(dims) == 3
            splits = 4;
        else
            splits = 7;
        end
    end
    % 'PlotMethod' given?
    [NVpairs,plotMethod] = readNameValuePair(NVpairs,'PlotMethod',{},'outer');

    % check name-value pairs
    if CHECKS_ENABLED
        if ~isnumeric(splits) || ~isscalar(splits)
            throw(CORAerror("CORA:wrongValue","name-value pair 'Splits'",'numeric scalar'))
        end
        if ~ismember(plotMethod,{'inner','outer'})
            throw(CORAerror("CORA:wrongValue","name-value pair 'PlotMethod'","'outer','inner'"))
        end
    end
end

function [ls,dims] = aux_preprocess(ls,dims)
    % preprocess

    if length(dims) == 1
        % add artificial dimension at 2nd dimension
        dim_old = 1:length(ls.vars);
        dim_old = 2 + dim_old; % shift
        dim_old(dims) = 1;
        ls = lift_(ls, length(dim_old)+2, dim_old);
        dims = [1;2];
    end
end

function han = aux_plotNd(ls,dims,NVpairs,splits,plotMethod)
    % different types of level sets
    if strcmp(ls.compOp,'==')
    
        % different methods for the different dimensions
        if length(dims) == 2
            han = aux_plot2Dcontour(ls,dims,NVpairs);
        else
            [res,ind] = aux_isSolvable(ls,dims);
    
            if res
                han = aux_plot3Dsolvable(ls,dims,ind,NVpairs);
            else
                han = aux_plot3Dsplit(ls,dims,splits,plotMethod,NVpairs); 
            end
        end
    
    else
        % different methods for the differnt dimensions
        if length(dims) == 2
            han = aux_plot2Dsplit(ls,dims,splits,plotMethod,NVpairs);
        else
            han = aux_plot3Dsplit(ls,dims,splits,plotMethod,NVpairs); 
        end 
    end
end

function han = aux_plot2Dcontour(obj,dims,type)
% plot 2D level set using Matlabs contour plot function

    % re-read plotOptions, since always plot called
    type = readPlotOptions(type,'contour');

    % get limits of current plot
    [xLim,yLim,~] = getUnboundedAxisLimits();
    
    % substitute all remaining entries with zero
    p = zeros(obj.dim,1);

    % generate contour plot
    N = 30;
    x = xLim(1):(xLim(2)-xLim(1))/N:xLim(2);
    y = yLim(1):(yLim(2)-yLim(1))/N:yLim(2);

    [X,Y] = meshgrid(x,y);
    Z = zeros(size(X));

    for i = 1:size(Z,1)
        for j = 1:size(Z,2)
            p_ = p;
            p_(dims) = [X(i,j);Y(i,j)];
            Z(i,j) = obj.funHan(p_);
        end
    end

    % level at which contour is plotted: always at z = 0
    level = [0 0];
    [~,han] = contour(X,Y,Z,level,type{:});
    updateColorIndex; % does not get updated for contour plots
end

function han = aux_plot2Dsplit(obj,dims,splits,plotMethod,NVpairs)
% plot 2D level set by splitting the plot area 

    % re-read plotOptions, since always fill called
    NVpairs = readPlotOptions(NVpairs,'fill');


    % get limits of current plot
    [xLim,yLim,~] = getUnboundedAxisLimits();
    space = interval([xLim(1);yLim(1)],[xLim(2);yLim(2)]);

    % determine subspaces to be plotted
    subSpaces = aux_refineSpace(obj,dims,space,splits,plotMethod);

    % plot ---

    % read plot settings (TODO: move to new function (also for reachset etc.)?)
    holdStatus = ishold;
    if ~holdStatus
        plot(NaN,NaN,'HandleVisibility','off');
        % reset color index (before readPlotOptions!)
        set(gca(),'ColorOrderIndex',1);

        % reset limits
        xlim(xLim);
        ylim(yLim);
    end
    hold on;
    ax = gca();
    oldColorIndex = ax.ColorOrderIndex;

    % plot sets (TODO: unify)
    for i=1:length(subSpaces)
        han_i = plot(subSpaces{i},[1,2],NVpairs{:});
        if i == 1
            han = han_i;
            NVpairs = [NVpairs, {'HandleVisibility','off'}];
        end
    end
    if isempty(subSpaces)
        [NVpairs,facecolor] = readNameValuePair(NVpairs,'FaceColor');
        han = fill(nan,nan,facecolor,NVpairs{:});
    end

    % restore plot settings
    updateColorIndex(oldColorIndex);
    if ~holdStatus
        hold off;
    end
end

function subSpaces = aux_refineSpace(obj,dims,space,splits,plotMethod)
    % refine space and return included interval subspaces
    subSpaces = {};

    % init n-dim space (with all other dimenions = 0)
    nSpace = interval(zeros(dim(obj),1));
    nSpace(dims) = space;

    % range bounding using fast interval arithmetic
    val = obj.funHan(nSpace);

    if all(val.sup <= 0) % compOp always '<=' or '<' here
        % if all equations satisfied
        subSpaces{1} = space;
    elseif any(val.inf > 0)
        % if any equation is unsatisfiable
        % don't include subspace
    else
        % not (yet) decidedable
        if splits == 0
            switch plotMethod
                case 'outer'
                    % plot
                    subSpaces{1} = space;
                case 'inner'
                    % don't plot
            end
        else
            % refine subspace
            partSpace = partition(space, 2);
            
            % check subspaces
            for i = 1:length(partSpace)
                subSpaces_i = aux_refineSpace(obj,dims,partSpace{i},splits-1,plotMethod);
                subSpaces = [subSpaces,subSpaces_i];
            end
        end
    end
end

function han = aux_plot3Dsolvable(obj,dims,ind,type)
% plot 3D level set by solving for one variable

    % re-read plotOptions, since always surf ~ fill called
    type = readPlotOptions(type,'fill');

    % get limits of current plot
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    lim = [xLim;yLim;zLim];
    
    % substitute all remaining entries with zero
    p = zeros(obj.dim,1);
    
    % get function handle
    f = obj.solved{ind}.funHan{1}.eq;
    
    % generate grid
    N = 100; 
    ind_ = setdiff(1:3,ind);
    
    dx = (lim(ind_(1),2)-lim(ind_(1),1))/N;
    dy = (lim(ind_(2),2)-lim(ind_(2),1))/N;
    
    dx_ = dx/2; dy_ = dy/2;
    
    x = lim(ind_(1),1)+dx_:dx:lim(ind_(1),2)-dx_;
    y = lim(ind_(2),1)+dy_:dy:lim(ind_(2),2)-dy_;
    
    [X,Y] = meshgrid(x,y);
    Z = zeros(size(X));
    
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            p_ = p;
            p_(dims(ind_)) = [X(i,j);Y(i,j)];
            Z(i,j) = f(p_);
        end
    end
    
    if ind == 1
        han = surf(Z,X,Y,type{:});
    elseif ind == 2
        han = surf(X,Z,Y,type{:});
    else
        han = surf(X,Y,Z,type{:});
    end

end

function han = aux_plot3Dsplit(obj,dims,splits,plotMethod,NVpairs)
% plot 3D level set by splitting the plot area

    % re-read plotOptions, since always fill called
    NVpairs = readPlotOptions(NVpairs,'fill');

    % get limits of current plot
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    space = interval([xLim(1);yLim(1);zLim(1)],[xLim(2);yLim(2);zLim(2)]);

    % determine subspaces to be plotted
    subSpaces = aux_refineSpace(obj,dims,space,splits,plotMethod);

    % plot ---

    % read plot settings (TODO: move to new function (also for reachset etc.)?)
    holdStatus = ishold;
    if ~holdStatus
        plot(NaN,NaN,'HandleVisibility','off');
        % reset color index (before readPlotOptions!)
        set(gca(),'ColorOrderIndex',1);

        % reset limits
        xlim(xLim);
        ylim(yLim);
        zlim(zLim);
    end
    hold on;
    ax = gca();
    oldColorIndex = ax.ColorOrderIndex;

    % plot sets
    for i=1:length(subSpaces)
        han_i = plot(subSpaces{i},[1,2,3],NVpairs{:});
        if i == 1
            han = han_i;
            NVpairs = [NVpairs, {'HandleVisibility','off'}];
        end
    end
    if isempty(subSpaces)
        [NVpairs,facecolor] = readNameValuePair(NVpairs,'FaceColor');
        han = fill3(nan,nan,nan,facecolor,NVpairs{:});
    end

    % restore plot settings
    updateColorIndex(oldColorIndex);
    if ~holdStatus
        hold off;
    end
end

function [res,ind] = aux_isSolvable(obj,dims)
% check if the level set equation is solvable for one variable

    res = false;
    
    if ~obj.solvable
       return; 
    end

    % find variable for which the equations are solved
    ind = 0;
    
    for i = 1:length(dims)
        if length(obj.solved{dims(i)}.eq) == 1
            ind = i;
            break;
        end
    end
    
    if ind ~= 0
       res = true;
    end
end

% ------------------------------ END OF CODE ------------------------------
