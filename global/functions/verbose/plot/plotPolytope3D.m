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
% See also: plotPolygon, plotPoints

% Authors:       Niklas Kochdumper
% Written:       02-December-2020
% Last update:   18-April-2024 (TL, fix color order and legend entries)
%                10-December-2024 (TL, output all graphic handles)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default
NVpairs = {'Color',nextcolor};

% read plot options if provided
if ~isempty(varargin)
    NVpairs = readPlotOptions(varargin);
end
% readout 'FaceColor' to decide plot/fill call where necessary
[~,facecolor] = readNameValuePair(NVpairs,'FaceColor');

% check if vertex array is not empty
if isempty(V)
    
    % plot dummy set...
    han = plot3(nan,nan,nan,NVpairs{:});

else
    % compute triangulation using the convex hull
    V = V';
    try
        [k,~] = convhulln(V);
    catch ME
        % not really 3D...

        try % try with plotPolygon            
            han = plotPolygon(V',NVpairs{:},'ConvHull',true);
            return
        catch 
            rethrow(ME)
        end
    end
    
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

    % save color order index
    ax = gca();
    oldColorOrderIndex = ax.ColorOrderIndex;

    hold on;
    
    % loop over all facets
    han = [];
    for i = 1:length(K)
        for j = 1:length(K{i})
            
            % compute vertex connection using the convex hull
            B = gramSchmidt(D{i}');
            vert = V(K{i}{j},:);
            vert_ = B'*vert';
            ind = convhull(vert_(2:end,:)');
        
            % plot the facet
            if isempty(facecolor) || strcmp(facecolor,'none')
                han_ij = plot3(vert(ind,1),vert(ind,2),vert(ind,3), NVpairs{:});
            else
                han_ij = fill3(vert(ind,1),vert(ind,2),vert(ind,3), facecolor, NVpairs{:});
            end
            han = [han;han_ij];

            if i == 1 && j == 1
                NVpairs = [NVpairs, {'HandleVisibility', 'off'}];
            end
        end
    end

    % reset hold status
    if ~holdStatus
        hold('off');
    end
    
    % update color order index
    updateColorIndex(oldColorOrderIndex);
    
    % show z-axis if currently not visible
    [az,~] = view();
    if az == 0
        view(-35, 30);
    end
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
