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

% project the object to the 2D-subspace
cZ = project(cZ,dims);

% plot modes: standard (1), template (2), splits (3)
if mode == 1
    han = plotStandard(cZ,plotOptions);
elseif mode == 2
    han = plotSplit(cZ,splits,plotOptions);
elseif mode == 3
    han = plotTemplate(cZ,numDir,plotOptions);
end


end


% Auxiliary Functions -----------------------------------------------------

function han = plotStandard(cZ,plotOptions)

% calculate vertices of the constrained zonotope
if isempty(cZ.A)
    V = polygon(zonotope(cZ.Z));
else
    V = vertices(cZ);
end

% plot a polygon that connects the vertices
han = plotPolygon(V,plotOptions{:});

end


function han = plotSplit(cZ,splits,plotOptions)

% recursively split the constrained zonotope
list = {cZ};

for i = 1:splits
   
    listTemp = cell(length(list)*2,1);
    counter = 1;
    
    % loop over all sets at the current recursion level
    for j = 1:length(list)
        
       % calculate radius of the interval over-approximation as a heuristic
       % indicating which dimension should be best splitted
       inter = interval(list{j});
       r = rad(inter);
        
       % split the set
       if r(1) > r(2)
           temp = split(list{j},1);
       else
           temp = split(list{j},2);
       end
       
       % update variables
       listTemp{counter} = temp{1};
       listTemp{counter+1} = temp{2};
       counter = counter + 2;
    end
    
    list = listTemp;
end

for i = 1:length(list)
    % over-approximate the splitted sets with intervals
    han = plot(interval(list{i}),[1,2],plotOptions{:}); hold on;
end


end


function han = plotTemplate(cZ,numDir,plotOptions)

% select directions for template polyhedron 
angles = linspace(0,360,numDir+1);
angles = angles(1:end-1);
angles = deg2rad(angles);

C = zeros(2,length(angles));

for i = 1:length(angles)
	C(:,i) = [cos(angles(i));sin(angles(i))];
end

% calculate the upper bounds along the directions
d = zeros(length(angles),1);

for i = 1:length(d)
    d(i) = supportFunc(cZ,C(:,i),'upper');
end

% construct template polyhedron
poly = mptPolytope(C',d);

% plot the template polyhedron
han = plot(poly,[1,2],plotOptions{:});

end

%------------- END OF CODE --------------