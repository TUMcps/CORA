function han = plot(cZ,varargin)
% plot - Plot a 2D-projection of a constrained zonotope object
%
% Syntax:  
%    plot(cZ)
%    plot(cZ,dims,type)
%
% Inputs:
%    cZ - constrained zonotope object
%    dims - (optional) dimensions of the projection
%    plotOptions - (optional) plot settings (LineSpec and name-value pairs)
%
% Outputs:
%    han - handle of figure
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    figure;
%    plot(cZ,[1,2],'r');
%
%    figure;
%    plot(cZ,[1,2],'b','Filled',true,'Splits',4);
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
% Last revision:---

%------------- BEGIN CODE --------------

% default settings
dims = [1,2];
plotOptions{1} = 'b';
mode = 1; % standard plot mode

% parse input arguments
if nargin >= 2 && ~isempty(varargin{1})
	dims = varargin{1}; 
end
if nargin >= 3 && ~isempty(varargin{2})
    plotOptions = varargin(2:end);
    % process linespec and Name-Value pairs
    [linespec,NVpairs] = readPlotOptions(plotOptions);
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar');
    [NVpairs,numDir] = readNameValuePair(NVpairs,'Template','isscalar');
    
    % error if both given
    if ~isempty(splits) && ~isempty(numDir)
        error("Choose either Splits or Template.");
    elseif ~isempty(splits)
        mode = 2;
    elseif ~isempty(numDir)
        mode = 3;
    end
    
    % unify with linespec
    plotOptions = [linespec,NVpairs];
end

% check dimension
if length(dims) < 2
    error('At least 2 dimensions have to be specified!');
elseif length(dims) > 3
    error('Only up to 3 dimensions can be plotted!');
end

% project the object to the 2D-subspace
cZ = project(cZ,dims);

% plot modes: standard (1), template (2), splits (3)
if mode == 1
    han = plotStandard(cZ,dims,plotOptions);
elseif mode == 2
    han = plotSplit(cZ,splits,dims,plotOptions);
elseif mode == 3
    han = plotTemplate(cZ,numDir,dims,plotOptions);
end

end


% Auxiliary Functions -----------------------------------------------------

function han = plotStandard(cZ,dims,plotOptions)

    % get dimensions after projection
    if length(dims) == 2
        dims = [1,2];
    else
        dims = [1,2,3];
    end

    if isempty(cZ.A)
        han = plot(zonotope(cZ.Z),dims,plotOptions{:});
    else
        
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
        if length(dims) == 2
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
        d(i) = supportFunc(cZ,C(:,i),'upper');
    end

    % construct template polyhedron
    poly = mptPolytope(C',d);

    % plot the template polyhedron
    han = plot(poly,dims_,plotOptions{:});

end

%------------- END OF CODE --------------