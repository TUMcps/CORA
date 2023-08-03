function han = plot(P,varargin)
% plot - plots a projection of a polytope
%
% Syntax:  
%    han = plot(P)
%    han = plot(P,dims)
%    han = plot(P,dims,type)
%
% Inputs:
%    P - mptPolytope object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Tobias Ladner
% Written:      19-October-2010
% Last update:  24-March-2015
%               17-March-2017
%               20-April-2018 (exception for empty sets)
%               25-May-2022 (TL: 1D Plotting)
%               22-November-2022 (TL: redundant VRep)
%               05-April-2023 (TL: clean up using plotPolygon)
%               26-July-2023 (TL: getUnboundedAxisLimits)
% Last revision:12-July-2023 (TL, restructure)

%------------- BEGIN CODE --------------

% 1. parse input
[P,dims,NVpairs] = aux_parseInput(P,varargin{:});

% 2. preprocess
[V,dims] = aux_preprocess(P,dims);

% 3. plot n-dimensional set
han = aux_plotNd(V,dims,NVpairs);

% 4. clear han
if nargout == 0
    clear han;
end

end

% Auxiliary functions -----------------------------------------------------

function [P,dims,NVpairs] = aux_parseInput(P,varargin)
    % parse input

    % default settings
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{P,'att','mptPolytope'};
                    {dims,'att','numeric',{'nonnan','vector','positive','integer'}}});
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % process linespec and Name-Value pairs
    NVpairs = readPlotOptions(varargin(2:end));
end

function [V,dims] = aux_preprocess(P,dims)
    % transform to 2D polytope if only 1D should be plotted
    if length(dims) == 1
        P = project(P,dims);
    end
     
    % check if polyhedron is bounded
    if ~isBounded(P.P)
        
        % get limits of current plot
        [xLim,yLim,zLim] = getUnboundedAxisLimits();
        
        % intersect with the current polytope
        if length(dims) == 2
            P = and_( ...
                project(P,dims), ...
                interval([xLim(1);yLim(1)],[xLim(2);yLim(2)]), ...
                'exact');
            dims = [1,2];
        else
            P = and_( ...
                project(P,dims), ...
                interval([xLim(1);yLim(1);zLim(1)], ...
                    [xLim(2);yLim(2),zLim(2)]), ...
                'exact');
            dims = [1,2,3];
        end
    end
        
    % compute vertices
    if ~P.P.irredundantVRep
        % solves issues with collinear points and convHull in plotPolygon
        P.P.computeHRep();
    end

    % project
    V = vertices(P);
    V = V(dims,:);
    dims = 1:size(V,1);
end

function han = aux_plotNd(V,dims,NVpairs)
    % plot projected vertices
    if isempty(V) 
        % plot empty set
        han = plotPolygon(zeros(length(dims), 0),NVpairs{:});

    elseif size(V,1) <= 2 % 1d, 2d
        han = plotPolygon(V(dims,:),NVpairs{:},'ConvHull',true);

    elseif size(V,1) == 3 % 3d
        han = plotPolytope3D(V(dims,:),NVpairs{:});
    end
end

%------------- END OF CODE --------------