function han = plotZono(cZ,varargin)
% plotZono - Visualizes a 2D-projection of the constraint zonotope and the
%            zonotope without constraints
%
% Syntax:  
%    plotZono(cZ)
%    plotZono(cZ,dims,plotOptZ,plotOptCon)
%
% Inputs:
%    cZ - constrained zonotope object
%    dims - (optional) dimensions of the projection
%    plotOptZ - (optional) cell-array containing the plot settings
%                 for the original zonotope
%    plotOptCon - (optional) cell-array containing the plot settings
%                     for the constrained zonotope
%
% Outputs:
%    han - handle of graphics object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    plotOptZ = {'r','LineWidth',2};
%    plotOptCon = {'b','EdgeColor','none'};
%    plotZono(cZ,[1,2],plotOptZ,plotOptCon);
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Author:       Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default settings
dims = [1,2];
plotOptZ = {'b'};
plotOptCon = {'r','EdgeColor','none'};

% parse input arguments
if nargin >= 2 && ~isempty(varargin{1})
	dims = varargin{1}; 
end
if nargin >= 3 && ~isempty(varargin{2})
	plotOptZ = varargin{2}; 
end
if nargin >= 4 && ~isempty(varargin{3})
	plotOptCon = varargin{3}; 
end

% plot the original zonotope
zono = zonotope(cZ.Z);
plot(zono,dims,plotOptZ{:});

% plot the constrained zonotope and return handle
hold on;
han = plot(cZ,dims,plotOptCon{:},'Filled',true);

%------------- END OF CODE --------------