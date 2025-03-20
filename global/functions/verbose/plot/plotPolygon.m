function han = plotPolygon(V,varargin)
% plotPolygon - plot a polygon defined by its vertices
%
% Syntax:
%    han = plotPolygon(V,varargin)
%
% Inputs:
%    V - matrix storing the polygon vertices
%    varargin - plot settings specified as linespec and name-value pairs
%    ('ConvHull', <true/false>) - whether convex hull of V should be computed
%    ('CloseRegions', <true/false>) - whether all given regions should closed
%    ('PlotBackground', <true/false>) - whether all given regions should closed by properly this function
%
% Outputs:
%    han - handle to the graphics object
%
% Example:
%    V = [1 0; 1 2; 0 3; -2 2; -3 0; 0 -1; 1 0]';
%    plotPolygon(V,'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plotPoints, plotPolytope3D

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       05-May-2020
% Last update:   15-July-2020 (MW, merge with plotFilledPolygon)
%                05-April-2023 (TL, generalized function)
%                11-July-2023 (TL, bug fix sets with holes and FaceColor)
%                12-July-2023 (TL, cut off infinity values at axis limits)
%                29-February-2024 (TL, fix ColorOrderIndex in filled 3d)
%                05-April-2024 (TL, added option to plot in background)
%                16-October-2024 (TL, fixes during contSet/plot restructuring)
%                17-December-2024 (TL, added 'CloseRegions')
%                12-February-2025 (TL, fixed plotting with multiple regions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default
NVpairs = {'Color',nextcolor};

% read plot options
if ~isempty(varargin)
    NVpairs = readPlotOptions(varargin);
end

% read additional name-value pairs
[NVpairs, doConvHull] = readNameValuePair(NVpairs, 'ConvHull', 'islogical', false);
[NVpairs, closeRegions] = readNameValuePair(NVpairs, 'CloseRegions', 'islogical', false);
[NVpairs, plotBackground] = readNameValuePair(NVpairs, 'PlotBackground', 'islogical', false);

% read positions
[NVpairs, V] = aux_positionAtXYZ(V, NVpairs);

% correct infinity values
V = aux_cutInfinityAtLimits(V);

% readout 'FaceColor' to decide plot/fill call where necessary
[NVpairs,facecolor] = readNameValuePair(NVpairs,'FaceColor');
hasFaceColor = ~isempty(facecolor) && ~strcmp(facecolor,'none');

% plot
if isempty(V)
    % plot empty (visible in legend)
    han = aux_plotEmpty(hasFaceColor, facecolor, NVpairs);

elseif size(V, 1) == 1
    % should be 2-dimensional after aux_positionAtXYZ() call
    throw(CORAerror("CORA:plotProperties", 2))

elseif size(V, 2) == 1
    % plot single point
    han = aux_plotSinglePoint(V,NVpairs,hasFaceColor,facecolor);

elseif size(V, 1) == 2
    % plot polygon in 2d
    han = aux_plot2D(V,NVpairs,doConvHull,hasFaceColor,facecolor,closeRegions);

elseif size(V, 1) == 3
    % plot polygon in 3d
    han = aux_plot3D(V,NVpairs,doConvHull,hasFaceColor,facecolor);

end

if plotBackground
    uistack(han, 'bottom');
end

if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [NVpairs, V] = aux_positionAtXYZ(V, NVpairs)
    
    % read XPos, YPos, ZPos
    [NVpairs,xpos] = readNameValuePair(NVpairs,'XPos',{'isnumeric','isscalar'});
    [NVpairs,ypos] = readNameValuePair(NVpairs,'YPos',{'isnumeric','isscalar'});
    [NVpairs,zpos] = readNameValuePair(NVpairs,'ZPos',{'isnumeric','isscalar'});

    % set default values
    [dims, n] = size(V);
    if dims == 1 && ~isempty(zpos)
        % if dim=1 and z is given, then specify either xpos or ypos
        if isempty(xpos) && isempty(ypos)
            ypos = 0;
        end
    end
    if dims == 1 && isempty(xpos) && isempty(ypos)
        ypos = 0;
    end

    % add respective dimensions to V
    if ~isempty(xpos)
        V = [xpos * ones(1, n); V];
    end
    if ~isempty(ypos)
        V = [V(1, :); ypos * ones(1, n); V(2:end, :)];
    end
    if ~isempty(zpos)
        V = [V; zpos * ones(1, n)];
    end

    % check dimensions to plot
    [dims, ~] = size(V);
    if dims < 2
        % this should not happen, just in case
        throw(CORAerror('CORA:specialError', 'Not enough dimensions specified.'))
    end
    if dims > 3
        throw(CORAerror('CORA:specialError', "Too many dimensions to plot. " + ...
            "Given points V along with optional parameters 'XPos', 'YPos', and 'ZPos' must not exceed 3!"))
    end 
end

function V = aux_cutInfinityAtLimits(V)
    % infinity values are cut off at axis limits
    % or smallest/largest value of V of respective dimension,
    % whichever is further 'outside'

    if any(isinf(V),'all')
        [xLim,yLim,zLim] = getUnboundedAxisLimits(V);

        % x-axis
        V(1,V(1,:) == -inf) = min([V(1,V(1,:) ~= -inf), xLim(1)]);
        V(1,V(1,:) == inf) = max([V(1,V(1,:) ~= inf), xLim(2)]);

        % y-axis
        V(2,V(2,:) == -inf) = min([V(2,V(2,:) ~= -inf), yLim(1)]);
        V(2,V(2,:) == inf) = max([V(2,V(2,:) ~= inf), yLim(2)]);
    
        if size(V,1) == 3
            % z-axis
            V(3,V(3,:) == -inf) = min([V(3,V(3,:) ~= -inf), zLim(1)]);
            V(3,V(3,:) == inf) = max([V(3,V(3,:) ~= inf), zLim(2)]);
        end
    end
end

function aux_show3dAxis()
    % show z-axis if currently not visible
    [az,~] = view();
    if az == 0
        view(-35, 30);
    end
end

function han = aux_plotEmpty(hasFaceColor, facecolor, NVpairs)
    if hasFaceColor
        han = fill(nan, nan, facecolor, NVpairs{:});
    else
        han = plot(nan, nan, NVpairs{:});
    end
end

function V = aux_convHull(V)
    % compute convex hull
    try
        ind = convhull(V(1,:),V(2,:));
        V = V(:,ind);
    catch
        if rank(V,1) <= 2
            % plot line, convhull fails if only collinear points are specified
        else
            throw(CORAerror('CORA:specialError','Plotting the set failed while constructing the convex hull.'));
        end
    end
end

function han = aux_plotSinglePoint(V,NVpairs,hasFaceColor,facecolor)
    % plot single point

    % check hold status and reset if required
    holdStatus = ishold;
    
    % plot empty (visible in legend)
    han = aux_plotEmpty(hasFaceColor, facecolor, NVpairs);
    % make scatter plot itself invisible
    NVpairs = [NVpairs, {'HandleVisibility', 'off'}];
    
    % correct color in NVpairs for scatter plot
    
    % EdgeColor > Color
    [NVpairs,color] = readNameValuePair(NVpairs, 'Color');
    [NVpairs,edgecolor] = readNameValuePair(NVpairs,'EdgeColor');
    NVpairs{end+1} = 'MarkerEdgeColor';
    if ~isempty(edgecolor)
        NVpairs{end+1} = edgecolor;
    else
        NVpairs{end+1} = color;
    end
    
    % FaceColor and marker
    if hasFaceColor
        % filled point
        NVpairs = [NVpairs,{'Marker','.','MarkerFaceColor',facecolor}];
    else
        % open circle
        NVpairs = [NVpairs,{'Marker','o'}];
    end

    % remove line style
    NVpairs = [{'LineStyle','none'},NVpairs];

    % set 'hold on' for actual plot to also show legend entry above
    hold on;
    
    if size(V, 1) == 2 
        % plot point in 2d
        plot(V(1,:), V(2,:),NVpairs{:});
    elseif size(V, 1) == 3 
        % plot point in 3d
        plot3(V(1,:),V(2,:),V(3,:),NVpairs{:});
        aux_show3dAxis()
    end

    % reset hold status
    if holdStatus
        hold on;
    else
        hold off;
    end
end

function han = aux_plot2D(V,NVpairs,doConvHull,hasFaceColor,facecolor,closeRegions)
    % plot in 2D
    
    if doConvHull
        % compute convex hull
        V = aux_convHull(V);
    end
    
    % axis are reset even if they were set manually; correct afterwards
    xmode = xlim("mode");
    if xmode == "manual"
        x = xlim;
    end
    ymode = ylim("mode");
    if ymode == "manual"
        y = ylim;
    end
    
    if ~hasFaceColor
        if closeRegions
            % close each region
            V = aux_closeRegions(V);
        end

        % make line plot
        han = plot(V(1,:), V(2,:), NVpairs{:});
    
    else
        % make filled plot

        % make sure it is a single region
        pgon = polygon(V);
        if pgon.set.NumRegions > 1
            % plot each individually
            han = plotMultipleSetsAsOne( ...
                arrayfun(@(pshape) polygon(pshape), pgon.set.regions,'UniformOutput',false), ...
                1:2, ... % dims
                [NVpairs(:)',{'FaceColor',facecolor,'CloseRegions',closeRegions}]);
            return;
        end

        % plot single region
        V = aux_sortMultipleRegionsAndHoles(V);
    
        % for the set with multiple regions/holes to fill correctly,
        % we need to remove the nan values
        idxNan = any(isnan(V),1);
        if nnz(idxNan) > 0
            % store start point of each region
            regStart = V(:,idxNan(2:end));
    
            % remove nan value
            V(:,any(isnan(V))) = [];
    
            % jump back to starting point V(:,1) to avoid line fragments
            V = [V,fliplr(regStart)];            
        end
    
        % sometimes the color index does not get increased automatically
        ax = gca();
        cidx = ax.ColorOrderIndex;
    
        % plot with facecolor
        han = fill(V(1,:), V(2,:), facecolor, NVpairs{:});
    
        if cidx == ax.ColorOrderIndex
            % update color index if it hasn't changed
            updateColorIndex;
        end
    end
    
    % reset axis mode
    if xmode == "manual"
        xlim(x);
    end
    if ymode == "manual"
        ylim(y);
    end
end

function V_closed = aux_closeRegions(V)

    % find nan values
    nanIdx = find(any(isnan(V),1));

    % find nan values
    regIdx = [0,nanIdx,size(V,2)+1];
    
    V_closed = [];
    for i=1:(numel(regIdx)-1)
        % read region
        V_reg = V(:,(regIdx(i)+1):(regIdx(i+1)-1));

        % check if closed
        if size(V_reg,2) > 1 && ~all(withinTol(V_reg(:,1),V_reg(:,2)),"all")
            % close region
            V_reg = [V_reg,V_reg(:,1)];
        end

        % append to all regions
        V_closed = [V_closed nan(2,1) V_reg];
    end

end

function V = aux_sortMultipleRegionsAndHoles(V)

% find region/hole separators (nan)
idxNan = find(any(isnan(V),1));
if isempty(idxNan)
    % single region, no hole
    return
end

% remove last point (= first point)
if all(withinTol(V(:,1),V(:,end)))
    V(:,end) = [];
end

% split vertices into regions
idxNan = [0,idxNan,size(V,2)+1];
Vs = arrayfun(@(i) V(:,(idxNan(i)+1):(idxNan(i+1)-1)), 1:(numel(idxNan)-1),'UniformOutput',false);

% determine largest region
[~,idx] = max(cellfun(@(V) vecnorm(max(V)-min(V)), Vs));
V_main = Vs{idx};
n = size(V,1);

% remove main region from regions
Vs(idx) = [];

% inserted regions/holes at smallest distance to boundary of main region
for i = 1:numel(Vs)
    % extract
    V_i = Vs{i};

    % compute minimum of each vertex with other vertex in main region
    V_i_3d = reshape(V_i, n, 1, []);
    diff = sum(abs(V_main-V_i_3d),1);
    [~, idxV] = min(min(diff, [], 3), [], 2);
    [~, idxH] = min(min(diff, [], 2), [], 3);

    V_main = [V_main(:, 1:idxV) ...
        V_i(:, idxH:end), V_i(:, 1:idxH), ...
        V_main(:, idxV:end)];
end


% obtain final vertices
V = V_main;

end


function han = aux_plot3D(V,NVpairs,doConvHull,hasFaceColor,facecolor)
% plot in 3D

if doConvHull
    % compute convex hull
    V = aux_convHull(V);
end

if ~hasFaceColor
    han = plot3(V(1,:),V(2,:),V(3,:),NVpairs{:}); 
else
    % sometimes the color index does not get increased automatically
    ax = gca();
    cidx = ax.ColorOrderIndex;

    han = fill3(V(1,:),V(2,:),V(3,:),facecolor,NVpairs{:}); 

    if cidx == ax.ColorOrderIndex
        % update color index if it hasn't changed
        updateColorIndex;
    end
end
aux_show3dAxis()

end

% ------------------------------ END OF CODE ------------------------------
