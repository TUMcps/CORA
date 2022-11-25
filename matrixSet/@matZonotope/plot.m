function plot(varargin)
% plot - Plots 2-dimensional projection of a matrix zonotope
%
% Syntax:  
%    plot(obj,dimensions)
%
% Inputs:
%    obj - matrix zonotope
%    dimensions - dimensions that should be projected (optional) 
%    linespec - plot style (optional) 
%
% Outputs:
%    none
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      22-June-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%convert from matrix zonotope to zonotope
Z=zonotope(varargin{1});
    
%plot zonotope
if nargin==1
    plot(Z);
elseif nargin==2
    plot(Z,varargin{2});
else
    plot(Z,varargin{2},varargin{3});
end

%------------- END OF CODE --------------