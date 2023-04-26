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

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      11-May-2018
% Last update:  15-July-2020 (MW, merge with plotFilled|Template|Split)
%               25-May-2022 (TL: 1D Plotting)
%               16-December-2022 (MW, add iterative method for 2D plots)
%               05-April-2023 (TL: clean up using plotPolygon)
% Last revision:---

%------------- BEGIN CODE --------------

% default settings
dims = setDefaultValues({[1,2]},varargin);
% standard plot mode
mode = 1;

% process linespec and Name-Value pairs
NVpairs = readPlotOptions(varargin(2:end));
[NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar');
[NVpairs,numDir] = readNameValuePair(NVpairs,'Template','isscalar');

% error if both given
if ~isempty(splits) && ~isempty(numDir)
    throw(CORAerror('CORA:specialError','Choose either Splits or Template.'));
elseif ~isempty(splits)
    mode = 2;
elseif ~isempty(numDir)
    mode = 3;
end

% check dimension
if length(dims) < 1
    throw(CORAerror('CORA:plotProperties',1));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end

% project the object to the 2D-subspace
cZ = project(cZ,dims);
dims = 1:length(dims);

% plot modes: standard (1), template (2), splits (3)
if mode == 1
    han = plotStandard(cZ,dims,NVpairs);
elseif mode == 2
    han = plotSplit(cZ,splits,dims,NVpairs);
elseif mode == 3
    han = plotTemplate(cZ,numDir,dims,NVpairs);
end

if nargout == 0
    clear han
end

end


% Auxiliary Functions -----------------------------------------------------

function han = plotStandard(cZ,dims,plotOptions)

    if isempty(cZ.A) || ( ~any(any(cZ.A)) && ~any(cZ.b) )
        han = plot(zonotope(cZ.Z),dims,plotOptions{:});
    elseif length(dims) == 2
        % 2D projection can be computed efficiently using support functions
        
        % compute vertices in projected dimensions
        V = projVertices(cZ,dims);

        if size(V,2) == 2
            % just a line... (does not work well with polygon/polyshape
            % class) -> instantiate zonotope and plot it
            c = 0.5*(V(:,1) + V(:,2)); 
            G = 0.5*(V(:,2) - V(:,1));
            han = plot(zonotope(c,G),dims,plotOptions{:});
        else
            % init polygon for plotting (vertices are already ordered
            % correctly)
            han = plotPolygon(V, plotOptions{:},'ConvHull',true);
        end

    else
        % other projections
        
        % convert to polytope
        poly = mptPolytope(cZ);
        
        % plot the polytope
        han = plot(poly,dims,plotOptions{:});        
    end
end

function han = plotSplit(cZ,splits,dims,plotOptions)

    % recursively split the constrained zonotope
    list = {cZ};

    for i = 1:splits

        listTemp = cell(length(list)*2,1);
        counter = 1;

        % loop over all sets at the current recursion level
        for j = 1:length(list)

           % calculate radius of the interval over-approximation as a
           % heuristic indicating which dimension should be best splitted
           inter = interval(list{j});
           r = rad(inter);

           % split the set
           [~,ind] = max(r);
           temp = split(list{j},ind);

           % update variables
           listTemp{counter} = temp{1};
           listTemp{counter+1} = temp{2};
           counter = counter + 2;
        end

        list = listTemp;
    end

    hold on;

    for i = 1:length(list)
        % over-approximate the splitted sets with intervals
        if length(dims) == 1
             han = plot(interval(list{i}),1,plotOptions{:});
        elseif length(dims) == 2
            han = plot(interval(list{i}),[1,2],plotOptions{:});
        else
            han = plot(interval(list{i}),[1,2,3],plotOptions{:});
        end
    end

end


function han = plotTemplate(cZ,numDir,dims,plotOptions)

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
        han = plotPolygon(V, plotOptions{:},'ConvHull',true);
    else

        % construct template polyhedron
        poly = polytope(C',d);
    
        % plot the template polyhedron
        han = plot(poly,dims_,plotOptions{:});

    end

end

%------------- END OF CODE --------------