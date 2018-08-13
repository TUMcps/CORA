function res = mptPolytope(Z)
% mptPolytope - Converts a zonotope object to a mptPolytope object
%
% Syntax:  
%    res = mptPolytope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    res - mptPolytope object
%
% Example: 
%    zono = zonotope(rand(2,5));
%    poly = mptPolytope(zono);
%
%    plot(poly,[1,2],'r');
%    figure
%    plot(zono,[1,2],'b');
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  vertices

% Author:       Niklas Kochdumper
% Written:      06-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = polytope(Z,'mpt');
    
%------------- END OF CODE --------------