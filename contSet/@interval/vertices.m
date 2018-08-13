function V = vertices(obj)
% vertices - Computes vertices of an interval object
%
% Syntax:  
%     V = vertices(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    V = vertices(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      24-Juli-2006 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%convert to zonotope 
Z = zonotope(obj);

%obtain vertices
V = vertices(Z);

%------------- END OF CODE --------------