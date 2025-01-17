function han = plot(S,varargin)
% plot - plots a projection of a contSet
%
% Syntax:
%    han = plot(S)
%    han = plot(S,dims)
%    han = plot(S,dims,varargin)
%
% Inputs:
%    S - contSet object
%    dims - (optional) dimensions for projection
%    varargin - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Authors:       Tobias Ladner
% Written:       14-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[S,dims,NVpairs,NVpairsInterval,NVpairsPolygon,NVpairsVertices] = aux_parseInput(S,varargin{:});

% process
[S,NVpairs] = aux_process(S,dims,NVpairs);

% call subfunction depending on number of dimensions to plot
n = dim(S);
if n == 1 % 1d
    han = plot1D(S,NVpairs,NVpairsInterval);
elseif n == 2 % 2d
    han = plot2D(S,NVpairs,NVpairsPolygon);
elseif n == 3 % 3d
    han = plot3D(S,NVpairs,NVpairsVertices);
else % unable to plot higher-dimensional sets
    throw(CORAerror('CORA:plotProperties',3))
end

% clear handle if not desired
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [S,dims,NVpairs,NVpairsInterval,NVpairsPolygon,NVpairsVertices] = aux_parseInput(S,varargin)
    % parse input

    % set default values
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{S,'att','contSet'};
                    {dims,'att','numeric',{'nonempty','integer','positive','vector'}}});
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % parse plot options
    NVpairs = varargin(2:end);

    % read additional parameters
    NVpairsInterval = {};
    NVpairsPolygon = {}; 
    NVpairsVertices = {};

    if isa(S,'polyZonotope')
        [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar',8);
        NVpairsInterval = {'split',splits};
        NVpairsPolygon = {splits};
        NVpairsVertices = {splits};
    elseif isa(S,'conPolyZono')
        [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar',8);
        NVpairsInterval = {'conZonotope'};
        NVpairsPolygon = {splits};
        NVpairsVertices = {splits};
    end
end

function [S,NVpairs] = aux_process(S,dims,NVpairs)
    % project set
    S = project(S,dims);

    % purpose for readPlotOptions
    purpose = 'none';

    % check if unbounded
    I = interval(S);
    if ~isBounded(I)
        % intersect with current plot axis
        S = aux_intersectWithAxisLimits(S,NVpairs);        

        % fill unbounded sets (fullspace, halfspace, ...)
        purpose = 'fill';
    end
    NVpairs = readPlotOptions(NVpairs,purpose);
end

function S = aux_intersectWithAxisLimits(S,NVpairs)
    % intersects unbounded sets with 

    % get quick estimate of vertices
    V = vertices_(interval(S));

    % consider given XPos, YPos, ZPos
    V = aux_positionVertices(V,NVpairs);
    n = size(V,1);

    % get size of current plot
    if n <= 2 % 1-2d
        [xLim,yLim] = getUnboundedAxisLimits(V);

        % init axis interval
        I_axis = interval([xLim(1);yLim(1)],[xLim(2);yLim(2)]);
        
    elseif n == 3 % 3d
        % get size of current plot
        [xLim,yLim,zLim] = getUnboundedAxisLimits(V);

        % init axis interval
        I_axis = interval([xLim(1);yLim(1);zLim(1)],[xLim(2);yLim(2);zLim(2)]);

    else % unable to plot higher-dimensional sets
        throw(CORAerror("CORA:plotProperties",3))
    end

    % project to given dimensions
    I_axis = project(I_axis,1:dim(S));

    % intersect with set
    S = and_(S,I_axis,'exact');

end

function V = aux_positionVertices(V,NVpairs)
    % read dimension
    [n,N] = size(V);
    z = zeros(1,N);

    % read XPos
    [NVpairs,xpos] = readNameValuePair(NVpairs,'XPos',{'isnumeric','isscalar'});
    if ~isempty(xpos)
        V = [z;V];
        n = n+1;
    end

    % read YPos
    [NVpairs,ypos] = readNameValuePair(NVpairs,'YPos',{'isnumeric','isscalar'});
    if ~isempty(ypos)
        V = [V(1,:);z,V(2:end)];
        n = n+1;
    end

    % read ZPos
    [~,zpos] = readNameValuePair(NVpairs,'ZPos',{'isnumeric','isscalar'});
    if ~isempty(zpos)
        if n == 1
            V = [V;z];
        end
        V = [V(1:2,:);z,V(3:end)];
    end
end

% ------------------------------ END OF CODE ------------------------------
