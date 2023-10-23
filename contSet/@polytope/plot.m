function han = plot(P,varargin)
% plot - Plots 2-dimensional projection of a polytope
%
% Syntax:
%    han = plot(P)
%    han = plot(P,dims)
%    han = plot(P,dims,type)
%
% Inputs:
%    P - polytope object
%    dims - (optional) dimensions of the projection
%    type - (optional) plot settings (LineSpec and name-value pairs)
%
% Outputs:
%    han - handle to the plotted graphics object
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       19-October-2010
% Last update:   24-March-2015
%                17-March-2017
%                20-April-2018 (exception for empty sets)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default settings
dims = setDefaultValues({[1,2]},varargin);

% check input arguments
inputArgsCheck({{P,'att','polytope'};
                {dims,'att','numeric',{'nonnan','vector','positive','integer'}}});

% process linespec and Name-Value pairs
NVpairs = readPlotOptions(varargin(2:end));

% check dimension
if length(dims) < 1
    throw(CORAerror('CORA:plotProperties',1));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end

% transform to 2D polytope if only 1D should be plotted
if length(dims) == 1
    P = project(P,dims);
    P = P.cartProd(0);
    dims = [1, 2];
end
 
% check if polyhedron is bounded
if ~isBounded(P)
    
    % get size of current plot
    xLim = get(gca,'Xlim');
    yLim = get(gca,'Ylim');
    
    % intersect with the current polytope
    if length(dims) == 2
        P = and_(project(P,dims),interval([xLim(1);yLim(1)], ...
                                           [xLim(2);yLim(2)]),'exact');
        dims = [1,2];
    else
        zLim = get(gca,'Zlim');
        P = and_(project(P,dims),interval([xLim(1);yLim(1);zLim(1)], ...
                                           [xLim(2);yLim(2);zLim(2)]),'exact');
        dims = [1,2,3];
    end
end
    
% compute vertices
% commented code below integrated to solve some issues
% % if ~P.irredundantVRep
% %     % solves issues with collinear points and convHull in plotPolygon
% %     P.computeHRep();
% % end
V = vertices(P);

% plot projected vertices
if length(dims) == 2
    han = plotPolygon(V(dims,:),NVpairs{:},'ConvHull',true);
elseif length(dims) == 3
    han = plotPolytope3D(V(dims,:),NVpairs{:});
end

if nargout == 0
    clear han;
end

% ------------------------------ END OF CODE ------------------------------
