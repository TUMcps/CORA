function plot(varargin)
% plot - Plots 2-dimensional projection of a pplPolytope
% Syntax:  
%    plot(P,dimensions,type)
%
% Inputs:
%    P - pplPolytope
%    dimensions - dimensions that should be projected (optional) 
%    type - plot type (optional) 
%
% Outputs:
%    none
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    P=varargin{1};
    dimensions=[1,2];
    type='frame';
    
%If two arguments are passed    
elseif nargin==2
    P=varargin{1};
    dimensions=varargin{2};
    type='frame';
    
%If too many arguments are passed
elseif nargin==3
    P=varargin{1};
    dimensions=varargin{2};   
    type=varargin{3};
end
 
%compute vertices
Pproj=project(P,dimensions);
V=vertices(Pproj);
%Plot convex hull of projected vertices
plot(V,type);

%------------- END OF CODE --------------