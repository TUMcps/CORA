function plotProj(varargin)
% plotProj - Plots 2-dimensional projections of a polytopes
%
% Syntax:  
%    plot(P,dimensions)
%
% Inputs:
%    P - polytope object
%    dimensions - dimensions that should be projected (optional) 
%
% Outputs:
%    none
%
% Example: ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 11-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    P=varargin{1};
    dimensions=[1,2];
    
%If two arguments are passed    
elseif nargin==2
    P=varargin{1};
    dimensions=varargin{2};
    
%If too many arguments are passed
else
    disp('Error: too many inputs');
    P=varargin{1};
    dimensions=varargin{2};    
end

%Compute projection
Pproj=projection(P,dimensions);

%convert to vertices object
V=extreme(Pproj)';
V=vertices(V);

%Plot vertices
plot(V);

%------------- END OF CODE --------------