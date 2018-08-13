function V = vertices(obj)
% vertices - computes the vertices of a pplPolytope
%
% Syntax:  
%    V = vertices(obj)
%
% Inputs:
%    obj - pplPolytope object
%
% Outputs:
%   V - vertices object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      19-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute vertices
V = ExtremePoints(-obj.C, obj.d);

%convert to vertices object
V = vertices(V');

%------------- END OF CODE --------------