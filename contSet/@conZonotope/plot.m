function plot(obj,varargin)
% plot - Plot a 2D-projection of a constrained zonotope object
%
% Syntax:  
%    plot(obj)
%    plot(obj,dim,color)
%
% Inputs:
%    obj - constrained zonotope object
%    dim - dimensions of the projection
%    plotOptions - plot settings specified as name-value pairs
%
% Outputs:
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    plot(cZono,[1,2],'r');
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
plotOptions = {'b'};

% parse input arguments
if nargin >= 2 && ~isempty(varargin{1})
   dim = varargin{1}; 
end
if nargin >= 3 && ~isempty(varargin{1})
   plotOptions = varargin(2:end); 
end

% project the object to the 2D-subspace
obj = project(obj,dim);

% calculate vertices of the constrained zonotope
v = vertices(obj);

% plot a polygon that connects the vertices
plot(v,plotOptions{:});

%------------- END OF CODE --------------