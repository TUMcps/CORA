function han = plotPolytope3D(V,varargin)
% plotPolytope3D - plot a polytope defined by its vertices in 3D
%
% Syntax:
%    han = plotPolytope3D(V,varargin)
%
% Inputs:
%    V - matrix storing the polytope vertices
%    varargin - plot settings specified as linespec and name-value pairs
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    zono = zonotope.generateRandom('Dimension',3,'NrGenerators',10);
%    V = vertices(zono);
%
%    plotPolytope3D(V,'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plotPolygon

% Authors:       Niklas Kochdumper
% Written:       02-December-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default
NVpairs = {'Color',colorblind('b')};

% read plot options if provided
if ~isempty(varargin)
    NVpairs = readPlotOptions(varargin);
end
% readout 'FaceColor' to decide plot/fill call where necessary
[~,facecolor] = readNameValuePair(NVpairs,'FaceColor');

% check if vertex array is not empty
if isempty(V)
    
    % plot dummy set...
    han = plot3(nan,nan,nan,'HandleVisibility','off');

else

    % compute triangulation using the convex hull
    V = V';
    [k,~] = convhulln(V);
    
    % compute normal vectors of facets
    C = zeros(size(k,1),3);
    
    for i = 1:size(k,1)
       C(i,:) = cross(V(k(i,1),:)-V(k(i,2),:),V(k(i,2),:)-V(k(i,3),:))';
       C(i,:) = C(i,:)/norm(C(i,:));
       if sum(C(i,:)) < 0
           C(i,:) = -C(i,:);
       end
    end
    
    % try to unify facets with identical normal vectors
    [C_,ind] = sortrows(C);
    k_ = k(ind,:);
    lower = 1;
    i = 1;
    D = {}; K = {};
    
    while lower <= size(C_,1)
        while i <= size(C_,1) && all(abs(C_(lower,:) - C_(i,:)) < 1e-14) 
            i = i + 1;
        end
        D{end+1} = C_(lower,:);
        K{end+1} = aux_getConnectedFacets(k_(lower:i-1,:));
        lower = i;
    end
    
    % set hold on
    holdStatus = ishold;
    if ~holdStatus
        % flush current axis
        plot(nan,nan,'HandleVisibility','off');
    end

    hold on;
    
    % loop over all facets
    for i = 1:length(K)
        for j = 1:length(K{i})
            
            % compute vertex connection using the convex hull
            B = gramSchmidt(D{i}');
            vert = V(K{i}{j},:);
            vert_ = B'*vert';
            ind = convhull(vert_(2:end,:)');
        
            % plot the facet
            if isempty(facecolor) || strcmp(facecolor,'none')
                han = plot3(vert(ind,1),vert(ind,2),vert(ind,3), NVpairs{:});
            else
                han = fill3(vert(ind,1),vert(ind,2),vert(ind,3), facecolor, NVpairs{:});
            end
        end
    end

    % reset hold status
    if ~holdStatus
        hold('off');
    end
    
    view([1,1,1]);
end
end


% Auxiliary functions -----------------------------------------------------

function vert = aux_getConnectedFacets(k)

    % initialize groups
    vert = mat2cell(k,ones(1,size(k,1)),3);
    groups = num2cell(1:1:size(k,1));
    
    cnt = 1;
    
    while cnt < length(groups)
        
        for i = cnt+1:length(groups)
            if any(ismember(vert{i},vert{cnt}))
                groups{cnt} = [groups{cnt};groups{i}];
                vert{cnt} = unique([vert{cnt},vert{i}]);
                groups{i} = [];
                vert{i} = [];
            end
        end
        
        cnt = cnt + 1;
        groups = groups(~cellfun('isempty',groups));
        vert = vert(~cellfun('isempty',vert));
    end
end

% ------------------------------ END OF CODE ------------------------------
