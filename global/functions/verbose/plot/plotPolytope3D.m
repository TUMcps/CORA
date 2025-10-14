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
%                14-July-2025 (TL, return in single graphics handle)
%                27-August-2025 (TL, bug fix, face color requires separate handles?)
%                11-September-2025 (TL, better handling of face color plots)
%                06-October-2025 (TL, bugfix degenerate convex sets)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse name-value pairs
if isempty(varargin)
    NVpairs = {'Color',nextcolor};
else
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
            han = plotPolygon(V',NVpairs{:},'ConvHull',true,'NVPAIRS_VALIDATED',true);
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
    
    % loop over all facets ---

    % preallocate
    N = sum(cellfun(@numel,K));
    % pre-allocate facets cell
    psCell = cell(N,1);
    cnt = 1;
    for i = 1:length(K)
        for j = 1:length(K{i})  
            % compute vertex connection using the convex hull
            B = gramSchmidt(D{i}');
            vert = V(K{i}{j},:);
            vert_ = B'*vert';
            ind = convhull(vert_(2:end,:)');

            % store respective vertices
            psCell{cnt} = [
                vert(ind,:);
                nan(1,3); % separating nan
            ];
            cnt = cnt+1; 
        end
    end

    % merge cell to matrix
    ps = cell2mat(psCell);

    % plot the facet
    if isempty(facecolor) || strcmp(facecolor,'none')
        han = plot3(ps(:,1),ps(:,2),ps(:,3), NVpairs{:});
    else
        % plot into one axis 
        % (doesn't show facecolor if there are nan values...)
        % -> plot each face individually using polytopes
        Ps = cellfun(@(ps) polytope(ps(1:end-1,:)'),psCell,'UniformOutput',false);
        han = plotMultipleSetsAsOne(Ps,1:3,[NVpairs,{'FaceColor',facecolor}]);
    end
    
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
