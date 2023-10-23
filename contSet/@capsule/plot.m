function han = plot(C, varargin)
% plot - plots a projection of a capsule
%
% Syntax:
%    han = plot(C)
%    han = plot(C,dims)
%    han = plot(C,dims,type)
%
% Inputs:
%    C - capsule object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example:
%    C = capsule([1; 1; 0], [0; 1; 1], 0.5);
%    plot(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Authors:       Matthias Althoff
% Written:       04-March-2019
% Last update:   25-May-2022 (TL, 1D Plotting)
%                05-April-2023 (TL, clean up using plotPolygon)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input arguments
[C, dims, NVpairs] = aux_parseInput(C, varargin{:});

% 2. plot n-dimensional set
han = aux_plotNd(C,dims,NVpairs);

% 3. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [C, dims, NVpairs] = aux_parseInput(C, varargin)
    % parse input arguments
    dims = setDefaultValues({[1, 2]}, varargin);
    % check input arguments
    inputArgsCheck({{C, 'att', 'capsule'}; ...
        {dims, 'att', 'numeric', {'nonempty', 'vector', 'integer', 'positive'}}});
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties', 1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties', 3));
    end
    
    % parse plot options
    NVpairs = readPlotOptions(varargin(2:end));
end

function han = aux_plotNd(C, dims, NVpairs)
    % plot n-dimensional set

    if length(dims) == 1 % 1d
        han = plot(interval(C), dims, NVpairs{:});
    
    elseif length(dims) == 2 % 2d
        han = aux_plot2d(C, dims, NVpairs);
    
    else % 3d
        han = aux_plot3d(C, dims, NVpairs);
    end
end

function han = aux_plot2d(C, dims, NVpairs)
    % plot 2-dimensional set

    % project capsule
    C = project(C, dims);
    
    % underapproximate capsule by polygon
    p = polygon(C);
    
    % plot and output the handle
    han = plotPolygon(p, NVpairs{:});
end
    
function han = aux_plot3d(C, dims, NVpairs)
    % plot 3-dimensional set

    % project capsule
    C = project(C, dims);

    if isempty(C.g)
        % generator is empty
    
        E = ellipsoid(eye(3)*C.r^2, C.c);
        Z = zonotope(E, 100, 'outer:norm');
        han = plot(Z, [1,2,3], NVpairs{:});
    
    else
        % generator is not empty
    
        % generate ellipsoids for spheres
        E1 = ellipsoid(eye(3)*C.r^2, C.c+C.g);
        E2 = ellipsoid(eye(3)*C.r^2, C.c-C.g);
    
        % enclose ellipsoid with a zonotope
        Z1 = zonotope(E1, 100, 'outer:norm');
        Z2 = zonotope(E2, 100, 'outer:norm');
    
        % compute vertices of the zonotopes
        V1 = vertices(Z1);
        V2 = vertices(Z2);
    
        % plot the convex hull
        han = plotPolytope3D([V1, V2], NVpairs{:});
    end
end

% ------------------------------ END OF CODE ------------------------------
