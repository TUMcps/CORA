function P = mptPolytope(cZ)
% mptPolytope - convert a constrained zonotope object to a mptPolytope
%    according to Prop. 3 in [1]
%
% Syntax:  
%    P = mptPolytope(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    P - mptPolytope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%    P = mptPolytope(cZ);
%
%    figure; hold on;
%    plot(cZ,[1,2],'FaceColor','r');
%    plot(P,[1,2],'g');
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
    
    % 2-dimensional set -> use efficient algorithm for vertex computation
    if dim(cZ) == 2
        V = vertices(cZ);
        P = mptPolytope(V');
        return;
    end

    % construct lifted zonotope 
    c = cZ.Z(:,1); 
    G = cZ.Z(:,2:end);
    
    Z = zonotope([c;-cZ.b],[G;cZ.A]);
    
    % compute halfspace reprsentation of lifted zonotope
    P = mptPolytope(Z);
    
    % extract halfspace representation for the constrained zonotope
    C = P.P.A;
    d = P.P.b;
    n = dim(cZ);
    
    % instantiate polytope
    P = mptPolytope(C(:,1:n),d);

%------------- END OF CODE --------------