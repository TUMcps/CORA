function res = mptPolytope(obj)
% mptPolytope - convert a constrained zonotope object to a mptPolytope
%               according to Prop. 3 in [1]
%
% Syntax:  
%    res = mptPolytope(obj)
%
% Inputs:
%    obj - conZonotope object
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
% See also: mptPolytope
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      13-May-2018
% Last update:  28-April-2019 (MA, code shortened)      
%               27-December-2020 (NK, use advanced method)
% Last revision:---

%------------- BEGIN CODE --------------
    
    % construct lifted zonotope 
    c = obj.Z(:,1); 
    G = obj.Z(:,2:end);
    
    zono = zonotope([c;-obj.b],[G;obj.A]);
    
    % compute halfspace reprsentation of lifted zonotope
    poly = mptPolytope(zono);
    
    % extract halfspace representation for the constrained zonotope
    C = poly.P.A;
    d = poly.P.b;
    n = dim(obj);
    
    res = mptPolytope(C(:,1:n),d);

%------------- END OF CODE --------------