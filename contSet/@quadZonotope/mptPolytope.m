function res = mptPolytope(qZ)
% mptPolytope - Computes a mptPolytope object that over-approximates the
%               quadZonotope object
%
% Syntax:  
%    res = mptPolytope(qZ)
%
% Inputs:
%    qZ - quadZonotope object
%
% Outputs:
%    res - mptPolytope object
%
% Example: 
%    qZ = quadZonotope([0;0],[1 0;0 1],[1 2;2 -1],[1;1],[0;0]);
%    poly = mptPolytope(qZ);
%
%    plot(zonotope(qZ),[1,2],'r');
%    figure
%    plot(poly,[1,2],'b');
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, interval

% Author:       Niklas Kochdumper
% Written:      06-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = enclosingPolytope(zonotope(qZ));

%------------- END OF CODE --------------