function han = plot(cZ,varargin)
% plot - plots a projection of a constrained zonotope
%
% Syntax:
%    han = plot(cZ)
%    han = plot(cZ,dims)
%    han = plot(cZ,dims,type)
%
% Inputs:
%    cZ - conZonotope object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%           additional Name-Value pairs:
%               <'Splits',splits> - number of splits for refinement
%               <'Template',dirs> - template directions
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    figure;
%    plot(cZ,[1,2],'r');
%
%    figure;
%    plot(cZ,[1,2],'FaceColor','b','Splits',4);
%
%    figure;
%    plot(cZ,[1,2],'r','Template',16);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       11-May-2018
% Last update:   15-July-2020 (MW, merge with plotFilled|Template|Split)
%                25-May-2022 (TL, 1D Plotting)
%                16-December-2022 (MW, add iterative method for 2D plots)
%                05-April-2023 (TL, clean up using plotPolygon)
%                27-April-2023 (VG, check if cZ is feasible)
%                09-May-2023 (TL, bugfix split plotting)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[cZ,dims,NVpairs,mode,splits,numDir] = aux_parseInput(cZ,varargin{:});

% 2. preprocess
[cZ,dims] = aux_preprocess(cZ,dims);

% 3. plot n-dimensional set
han = aux_plotNd(cZ,dims,NVpairs,mode,splits,numDir);

% 4. clear han
if nargout == 0
    clear han
end

end


% Auxiliary functions -----------------------------------------------------

function [cZ,dims,NVpairs,mode,splits,numDir] = aux_parseInput(cZ,varargin)
    % parse input

    % default settings
    dims = setDefaultValues({[1,2]},varargin);

    % check input args
    inputArgsCheck({{cZ, 'att', 'conZonotope'},
        {dims,'att','numeric',{'nonempty','vector','integer','positive'}}})

    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % process linespec and Name-Value pairs
    NVpairs = readPlotOptions(varargin(2:end));
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar');
    [NVpairs,numDir] = readNameValuePair(NVpairs,'Template','isscalar');
    
    % choose plot mode
    if ~isempty(splits) && ~isempty(numDir)
        % error if both given
        throw(CORAerror('CORA:specialError','Choose either Splits or Template.'));
    elseif ~isempty(splits)
        % plot mode 'split'
        mode = 2;
    elseif ~isempty(numDir)
        % plot mode 'template'
        mode = 3;
    else
        % default plot mode 'standard'
        mode = 1;
    end
end

function [cZ,dims] = aux_preprocess(cZ,dims)
    % project the object to the N=dimensional subspace
    cZ = project(cZ,dims);
    dims = 1:length(dims);
end

function han = aux_plotNd(cZ,dims,NVpairs,mode,splits,numDir)
    % check if constraints are feasible
    if representsa_(cZ,'emptySet',eps)
        % plot empty set
        han = plotPolygon(zeros(length(dims), 0), NVpairs{:});
    
    else        
        % plot modes: standard (1), template (2), splits (3)
        if mode == 2
            han = aux_plotSplit(cZ,splits,dims,NVpairs);
        elseif mode == 3
            han = aux_plotTemplate(cZ,numDir,dims,NVpairs);
        else
            % default plot mode
            han = aux_plotStandard(cZ,dims,NVpairs);
        end
    end
end

function han = aux_plotStandard(cZ,dims,NVpairs)

    if isempty(cZ.A) || ( ~any(any(cZ.A)) && ~any(cZ.b) )
        han = plot(zonotope(cZ.c,cZ.G),dims,NVpairs{:});
    elseif length(dims) == 2
        % 2D projection can be computed efficiently using support functions
        
        % compute vertices in projected dimensions
        V = projVertices(cZ,dims);

        if size(V,2) == 2
            % just a line... (does not work well with polygon/polyshape
            % class) -> instantiate zonotope and plot it
            c = 0.5*(V(:,1) + V(:,2)); 
            G = 0.5*(V(:,2) - V(:,1));
            han = plot(zonotope(c,G),dims,NVpairs{:});
        else
            % init polygon for plotting (vertices are already ordered
            % correctly)
            han = plotPolygon(V, NVpairs{:},'ConvHull',true);
        end

    else
        % other projections
        
        try 
            % convert to polytope
            poly = polytope(cZ);
        catch ME
            % if conversion to polytope fails (e.g., because the
            % constrained zonotope has too many generators), use vertices
            try
                V = vertices(cZ);
                % remove duplicates
                VV = uniquetol(V',1e-8,'ByRows',true)';
                han = plotPolytope3D(VV,NVpairs{:});
                return
            catch ME
                throw(CORAerror('CORA:specialIssue',...
                    ['Plotting of constrained zonotope failed! '...
                    'Neither polytope conversion nor vertex enumeration successful.']));
            end
        end
        
        % plot the polytope
        han = plot(poly,dims,NVpairs{:});
    end
end

function han = aux_plotSplit(cZ,splits,dims,NVpairs)

    % recursively split the constrained zonotope
    cZSplit = {cZ};

    for i = 1:splits
        % init
        listTemp = cell(length(cZSplit)*2,1);
        counter = 1;

        % loop over all sets at the current recursion level
        for j = 1:length(cZSplit)

           % calculate radius of the interval over-approximation as a
           % heuristic indicating which dimension should be best splitted
           inter = interval(cZSplit{j});
           r = rad(inter);

           % split the set
           [~,ind] = max(r);
           temp = split(cZSplit{j},ind);

           % update variables
           listTemp{counter} = temp{1};
           listTemp{counter+1} = temp{2};
           counter = counter + 2;
        end

        cZSplit = listTemp;
    end

    % convert splitted sets to intervals
    Is = cell(1,length(cZSplit));
    for i=1:length(cZSplit)
        Is{i} = interval(cZSplit{i});
    end

    % plot all sets as one
    han = plotMultipleSetsAsOne(Is,dims,NVpairs);
end


function han = aux_plotTemplate(cZ,numDir,dims,NVpairs)

    % select directions for template polyhedron
    if length(dims) == 2
        angles = linspace(0,360,numDir+1);
        angles = angles(1:end-1);
        angles = deg2rad(angles);

        C = zeros(2,length(angles));

        for i = 1:length(angles)
            C(:,i) = [cos(angles(i));sin(angles(i))];
        end

        dims_ = [1,2];
    else
        N = ceil(sqrt(numDir));
        theta = 2 * pi * linspace(0,1,N);
        phi = acos(1 - 2 * linspace(0,1,N));
        [phi,theta] = meshgrid(phi,theta);
        C = zeros(3,numel(phi));
        cnt = 1;

        for i = 1:size(phi,1)
            for j = 1:size(phi,2)
                C(1,cnt) = sin(phi(i,j)) .* cos(theta(i,j));
                C(2,cnt) = sin(phi(i,j)) .* sin(theta(i,j));
                C(3,cnt) = cos(phi(i,j));
                cnt = cnt + 1;
            end
        end

        dims_ = [1,2,3];
    end

    % calculate the upper bounds along the directions
    d = zeros(size(C,2),1);
    for i = 1:length(d)
        d(i) = supportFunc_(cZ,C(:,i),'upper');
    end

    % compute intersection of neighboring constraints in 2D case
    if length(dims) == 2

        % loop over all pairs constraints to compute vertices
        V = zeros(2,numDir);
        for i=1:numDir
            if i==numDir
                % last constraint with first constraint
                V(:,i) = C(:,[i,1])' \ d([i,1]);
            else
                V(:,i) = C(:,[i,i+1])' \ d(i:i+1);
            end
        end
        % remove duplicates using relative tolerance
        [V,IA] = uniquetol(V',1e-3,'ByRows',true);
        % re-order
        [~,order] = mink(IA,length(IA));
        V = V(order,:)';

        % init polygon for plotting
        han = plotPolygon(V, NVpairs{:},'ConvHull',true);
    else

        % construct template polyhedron
        poly = polytope(C',d);
    
        % plot the template polyhedron
        han = plot(poly,dims_,NVpairs{:});

    end

end

% ------------------------------ END OF CODE ------------------------------
