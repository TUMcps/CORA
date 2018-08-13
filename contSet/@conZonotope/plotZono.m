function plotZono(obj,varargin)
% plotZono - Visualizes a 2D-projection of the constraint zonotope and the
%            zonotope without constraints
%
% Syntax:  
%    plotZono(obj)
%    plotZono(obj,dim,plotOptZ,plotOptCon)
%
% Inputs:
%    obj - constrained zonotope object
%    dim - dimensions of the projection
%    plotOptZ - cell-array containing the plot settings for the original
%               zonotpope
%    plotOptCon - cell-array containing the plot settings for the
%                 constrained zonotope
%
% Outputs:
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%
%    plotOptZ = {'r','LineWidth',2};
%    plotOptCon = {'b','EdgeColor','none'};
%    plotZono(cZono,[1,2],plotOptZ,plotOptCon);
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ploot

% Author:       Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default settings
dim = [1,2];
plotOptZ = {'b'};
plotOptCon = {'r','EdgeColor','none'};

% parse input arguments
if nargin >= 2 && ~isempty(varargin{1})
   dim = varargin{1}; 
end
if nargin >= 3 && ~isempty(varargin{2})
   plotOptZ = varargin{2}; 
end
if nargin >= 4 && ~isempty(varargin{3})
   plotOptCon = varargin{3}; 
end

% plot the orignal zonotope
zono = zonotope(obj.Z);
plot(zono,dim,plotOptZ{:});

% plot the constrained zonotope
hold on
plotFilled(obj,dim,plotOptCon{:});

%------------- END OF CODE --------------