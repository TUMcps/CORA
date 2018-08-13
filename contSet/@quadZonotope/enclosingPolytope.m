function [P] = enclosingPolytope(varargin)
% polytope - Converts a zonotope to a polytope representation
%
% Syntax:  
%    [P] = polytope(Z)
%
% Inputs:
%    Z - polynomial zonotope object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    ---
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalhull,  vertices

% Author:       Matthias Althoff
% Written:      04-October-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

P = enclosingPolytope(zonotope(varargin{1}), varargin{2:end});

%------------- END OF CODE --------------