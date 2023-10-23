function han = plot(Z,varargin)
% plot - plots a projection of a zonotope
%
% Syntax:
%    han = plot(Z)
%    han = plot(Z,dims)
%    han = plot(Z,dims,type)
%
% Inputs:
%    Z - zonotope object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs),
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    Z = zonotope([1;0],[0.4 -1 0.4; -1 0.3 0]);
%    plot(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Authors:       Matthias Althoff
% Written:       27-July-2016
% Last update:   14-July-2020 (merge with plotFilled)
%                25-May-2022 (TL, 1D Plotting)
%                05-April-2023 (TL, clean up using plotPolygon)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[Z,dims,NVpairs] = aux_parseInput(Z,varargin{:});

% 2. preprocess
[Z,dims] = aux_preprocess(Z,dims);

% 3. plot n-dimensional set
han = aux_plotNd(Z,dims,NVpairs);

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [Z,dims,NVpairs] = aux_parseInput(Z,varargin)
    % parse input

    % set default values
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{Z,'att','zonotope'};
                    {dims,'att','numeric',{'nonempty','integer','positive','vector'}}});
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % parse plot options
    NVpairs = readPlotOptions(varargin(2:end));
end

function [Z,dims] = aux_preprocess(Z,dims)
    % preprocess

    % project zonotope
    Z = project(Z,dims);
    dims = 1:dim(Z);
end

function han = aux_plotNd(Z,dims,NVpairs)
    % plot n-dimensional set

    if length(dims) == 1 % 1d
        V = vertices(Z);
        han = plotPolygon(V,NVpairs{:});
    
    elseif length(dims) == 2 % 2d
        % convert zonotope to polygon
        p = polygon(Z);
    
        % plot and output the handle
        han = plotPolygon(p, NVpairs{:});
        
    else % 3d
        % compute vertices
        V = vertices(Z);
        
        % generate 3D plot
        if isempty(V)
            % plot empty set
            han = plotPolygon(zeros(length(dims), 0),NVpairs{:});
        else
            han = plotPolytope3D(V(dims,:),NVpairs{:});
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
