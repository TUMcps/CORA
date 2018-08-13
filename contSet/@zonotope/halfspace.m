function [obj] = halfspace(obj)
% halfspace - Generates halfspace representation of the zonotope
%
% Syntax:  
%    [obj] = halfspace(obj)
%
% Inputs:
%    obj - zonotope object
%
% Outputs:
%    obj - zonotope object
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
% Written:      07-May-2007 
% Last update:  06-April-2017
% Last revision:---

%------------- BEGIN CODE --------------

%convert zonotope to polytope and retrieve halfspace representation
P = polytope(obj);
H = get(P,'H');
K = get(P,'K');

%write to object structure
obj.halfspace.H=H;
obj.halfspace.K=K;
obj.halfspace.equations=length(K);

%------------- END OF CODE --------------