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

% Author: Matthias Althoff
% Written: 07-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%convert zonotope to polytope and retrieve halfspace representation
[H,K]=double(polytope(obj));
%write to object structure
obj.halfspace.H=H;
obj.halfspace.K=K;
obj.halfspace.equations=length(K);

%------------- END OF CODE --------------