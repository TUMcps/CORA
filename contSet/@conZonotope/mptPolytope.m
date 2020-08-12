function res = mptPolytope(obj)
% mptPolytope - convert a constrained zonotope object to a mptPolytope
%
% Syntax:  
%    res = mptPolytope(obj)
%
% Inputs:
%    obj - c-zonotope object
%
% Outputs:
%    res - mptPolytope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    poly = mptPolytope(cZono);
%
%    hold on
%    plot(cZono,[1,2],'r','Filled',true,'EdgeColor','none')
%    plot(poly,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      13-May-2018
% Last update:  28-April-2019 (MA, code shortened)
% Last revision:---

%------------- BEGIN CODE --------------
    
% extract vertices of the constrained zonotope
V = vertices(obj);

% construct MPTpolytope
res = mptPolytope(V');


%------------- END OF CODE --------------