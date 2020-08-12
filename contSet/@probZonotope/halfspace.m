function probZ = halfspace(probZ)
% halfspace - Generates halfspace representation of the probabilistic zonotope
%
% Syntax:  
%    probZ = halfspace(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    probZ - probabilistic zonotope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:        Matthias Althoff
% Written:       07-May-2007 
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

%convert zonotope to polytope and retrieve halfspace representation
[H,K]=double(polytope(probZ));
%write to object structure
probZ.halfspace.H=H;
probZ.halfspace.K=K;
probZ.halfspace.equations=length(K);

%------------- END OF CODE --------------