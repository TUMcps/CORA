function han = plotPolygon(V,varargin)
% plotPolygon - plot a polygon defined by its vertices
%
% Syntax:
%    han = plotPolygon(V,varargin)
%
% Inputs:
%    V - matrix storing the polygon vertices
%    varargin - plot settings specified as linespec and name-value pairs
%    ('ConvHull', <true/false>) - whether convex hull of V should be
%        computed
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
% See also:

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       05-May-2020
% Last update:   15-July-2020 (MW, merge with plotFilledPolygon)
%                05-April-2023 (TL, generalized function)
%                11-July-2023 (TL, bug fix sets with holes and FaceColor)
%                12-July-2023 (TL, cut off infinity values at axis limits)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default
NVpairs = {'Color',CORAcolor('CORA:next')};

% read plot options
if ~isempty(varargin)
    NVpairs = readPlotOptions(varargin);
end

% read doConvHull
[NVpairs, doConvHull] = readNameValuePair(NVpairs, 'ConvHull', 'islogical', false);

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
    % plot point

    % plot empty (visible in legend)
    han = aux_plotEmpty(hasFaceColor, facecolor, NVpairs);
    % make scatter plot itself invisible
    NVpairs{end+1} = 'HandleVisibility';
    NVpairs{end+1} = 'off';

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

    if size(V, 1) == 2 
        % plot point in 2d
        scatter(V(1,:), V(2,:),NVpairs{:});
    elseif size(V, 1) == 3 
        % plot point in 3d
        scatter3(V(1,:),V(2,:),V(3,:),NVpairs{:});
        aux_show3dAxis()
    end

elseif size(V, 1) == 2
    % plot polygon in 2d

    if doConvHull
        % compute convex hull
        V = aux_convHull(V);
    end

    % axis are reset even if they were set manually; correct afterwards
    xmode = xlim("mode");
    x = xlim;
    ymode = ylim("mode");
    y = ylim;

    if ~hasFaceColor
        % make line plot
        han = plot(V(1,:), V(2,:), NVpairs{:});

    else
        % make filled plot

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

elseif size(V, 1) == 3
    % plot polygon in 3d at given height

    if doConvHull
        % compute convex hull
        V = aux_convHull(V);
    end

    if ~hasFaceColor
        han = plot3(V(1,:),V(2,:),V(3,:),NVpairs{:}); 
    else
        han = fill3(V(1,:),V(2,:),V(3,:),facecolor,NVpairs{:}); 
    end
    aux_show3dAxis()

end

if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [NVpairs, V] = aux_positionAtXYZ(V, NVpairs)
    
    [NVpairs,xpos] = readNameValuePair(NVpairs,'XPos','isscalar');
    [NVpairs,ypos] = readNameValuePair(NVpairs,'YPos','isscalar');
    [NVpairs,zpos] = readNameValuePair(NVpairs,'ZPos','isscalar');

    % legacy 'Height'
    [NVpairs,height] = readNameValuePair(NVpairs,'Height','isscalar');
    if ~isempty(height)
        warning("CORA Warning: Plotting with 'Height' is deprecated. Use 'ZPos' instead.")
        if ~isempty(zpos)
            throw(CORAerror('CORA:notSupported', ...
                "Plotting with 'ZPos' and 'Height' specified is not allowed. Use 'ZPos' instead."))
        end
        zpos = height;
    end

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
        throw(CORAerror('CORA:specialError', "Too many dimensions to plot." + ...
            "Specified 'dims', 'XPos', 'YPos', and 'ZPos' must not exceed 3!"))
    end 
end

function V = aux_cutInfinityAtLimits(V)
    % infinity values are cut off at axis limits
    % or smallest/largest value of V of respective dimension,
    % whichever is further 'outside'

    [xLim,yLim,zLim] = getUnboundedAxisLimits(V);

    if any(isinf(V),'all')
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

% ------------------------------ END OF CODE ------------------------------
