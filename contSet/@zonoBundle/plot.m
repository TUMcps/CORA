function han = plot(zB,varargin)
% plot - plots an over-approximative projection of a zonotope bundle
%
% Syntax:
%    han = plot(zB)
%    han = plot(zB,dims)
%    han = plot(zB,dims,type)
%
% Inputs:
%    zB - zonoBundle object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    Z{1} = zonotope([0 1 1 -1;0 0 2 1]);
%    Z{2} = zonotope([2 1 1 -1;1 0 -1 -1]);
%    zB = zonoBundle(Z);
%
%    plot(zB,[1,2],'r','LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       09-November-2010 
% Last update:   13-February-2012
%                19-October-2015 (NK, accelerate by using polyshape class)
%                22-May-2022 (TL, 1D plots)
%                05-April-2023 (TL, clean up using plotPolygon)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[zB,dims,NVpairs] = aux_parseInput(zB,varargin{:});

% 2. preprocess
[zB,dims] = aux_preprocess(zB,dims);

% 3. plot n-dimensional set
han = aux_plotNd(zB,dims,NVpairs);

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [zB,dims,NVpairs] = aux_parseInput(zB,varargin)
    % parse input

    % default settings
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{zB,'att','zonoBundle'};
                    {dims,'att','numeric',{'nonempty','integer','positive','vector'}}});
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end

    % process plot options
    NVpairs = readPlotOptions(varargin(2:end));
end

function [zB,dims] = aux_preprocess(zB,dims)
    % preprocess

    if length(dims) == 1
        zB = project(zB,dims);
        % add zeros to 2nd dimension
        zB = cartProd_(zB,0,'exact');
        dims = [1;2];
    end
end

function han = aux_plotNd(zB,dims,NVpairs)
    % plot n-dimensional set
    
    if length(dims) == 2 % 1d, 2d
        han = aux_plot2d(zB,dims,NVpairs);   
        
    else % 3d
        han = aux_plot3d(zB,dims,NVpairs);
    end
end

function han = aux_plot2d(zB,dims,NVpairs)
    % plot 2-dimensional set

    % turn off warning
    w = warning;
    warning('off','all');

    % eps Z
    eps = 1e-7;
    Z_eps = zonotope([0;0], eps * eye(2));

    % compute polytopes
    P = cell(length(zB.parallelSets), 1);
    for i = 1:zB.parallelSets

        % delete zero generators
        Z = compact_(zB.Z{i},'zeros',eps);

        % project zonotope
        Z = project(Z,dims);

        % enlargen Z slightly:
        % intersections with polyshapes are inconsistent
        % if one polyshape is just a line. if we enlargen Z slightly, 
        % we can ignore these edge cases during intersection
        Z = Z + Z_eps;

        % convert to polyshape (Matlab built-in class)
        temp = polygon(Z);
        V = temp(:,2:end);
        P{i} = polyshape(V(1,:),V(2,:));
    end

    % intersect polytopes
    Pint = P{1};
    for i = 2:zB.parallelSets
        Pint = intersect(Pint,P{i});
    end

    % get vertices
    V = Pint.Vertices';

    % remove eps to get tight enclosure
    m = mean(V, 2);
    V0 = V - m;
    V0 = V0 - sign(V0) * eps;
    V = V0 + m;
    
    % test if all points are identical within tolerance
    inTol = all(withinTol(V(:,1),V,eps), "all");

    Pint = polyshape(V(1,:),V(2,:));
    if Pint.NumRegions > 0 && ~inTol
        % plot polytope
        han = plotPolygon(V,NVpairs{:},'ConvHull',true);
    else
        if size(V, 2) > 0 && inTol
            % just plot first point
            V = V(:, 1);
        end
        % plot line
        han = plotPolygon(V,NVpairs{:});
        
    end

    % reset warning state to previous setting
    warning(w);
end

function han = aux_plot3d(zB,dims,NVpairs)
    % plot 3-dimensional set
    
    % project to plotted dimensions
    zB = project(zB,dims);
    
    % compute vertices
    V = vertices(zB);
    
    % plot 3D polytope
    han = plotPolytope3D(V(dims,:),NVpairs{:}); 
end


% ------------------------------ END OF CODE ------------------------------
