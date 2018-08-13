function res = mptPolytope(Z)
% mptPolytope - Converts a zonotope object to a mptPolytope object
%
% Syntax:  
%    res = mptPolytope(Z)
%
% Inputs:
%    Z - zonotopeBundle object
%
% Outputs:
%    res - mptPolytope object
%
% Example: 
%    zono = zonotope(rand(2,5));
%    zB = zonotopeBundle({zono,[-1 0;0 1]*zono});
%    poly = mptPolytope(zB);
%
%    plot(poly,[1,2],'r');
%    figure
%    plot(zB,[1,2],'b');
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

    res = polytope(Z);
    
%------------- END OF CODE --------------