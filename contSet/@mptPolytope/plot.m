function plot(varargin)
% plot - Plots 2-dimensional projection of a mptPolytope
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
% Last update:  24-March-2015
%               17-March-2017
%               20-April-2018 (exception for empty sets)
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
 
%check if polyhedron is bounded; use own plotting capabilities for bounded
%polyhedra
if isBounded(P.P)
    %compute vertices
    V=vertices(P);
    %Plot convex hull of projected vertices
    if iscell(V)
        for i=1:length(V)
            if ~isempty(V{i})
                plot(vertices(V{i}(dimensions,:)),type);
            end
        end
    else
        if ~isempty(V)
            plot(vertices(V(dimensions,:)),type);
        end
    end
else
    %project polyhedron to dimensions for plotting
    dims = dimension(P);
    projDims = length(dimensions);
    
    %init projection matrix
    M = zeros(projDims,dims);
    
    for i = 1:projDims
        M(i,dimensions(i)) = 1;
    end
    Pproj = M*P;
    
    %plot projection
    plot(Pproj.P);
end


%------------- END OF CODE --------------