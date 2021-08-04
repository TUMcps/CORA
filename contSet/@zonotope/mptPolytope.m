function P = mptPolytope(Z)
% mptPolytope - converts a zonotope object to a mptPolytope object
%
% Syntax:  
%    P = mptPolytope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    P - mptPolytope object
%
% Example: 
%    zono = zonotope(rand(2,5));
%    poly = mptPolytope(zono);
%
%    figure; hold on;
%    plot(poly,[1,2],'r');
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

    P = polytope(Z,'mpt');
    
%------------- END OF CODE --------------