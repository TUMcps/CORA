function d = hausdorffDist(P1,P2)
% hausdorffDist - Calculates the Hausdorff distance between two polytopes
%                 or a polytope and a point
%
% Syntax:  
%    d = hausdorffDist(P1,P2)
%
% Inputs:
%    P1 - mptPolytope object
%    P2 - contSet object of single point
%
% Outputs:
%    d - Hausdorff distance
%
% Examples:
%    P1 = mptPolytope.generateRandom(2);
%    P2 = mptPolytope.generateRandom(2);
%   
%    d = hausdorffDist(P1,P2)
%
%    figure; hold on;
%    plot(P1,[1,2],'r');
%    plot(P2,[1,2],'b');
%
% References: 
%    [1] S. Koenig, "Computational Aspects of the Hausdorff Distancein 
%        Unbounded Dimension", 2018
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadZonotope/hausdorffDist

% Author:       Niklas Kochdumper
% Written:      05-September-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    d = 0;

    % get mptPolytope object
    if ~isa(P1,'mptPolytope')
       temp = P1; P1 = P2; P2 = temp;
    end
    
    % differnt cases for different types of sets
    if isnumeric(P2)
        
        for i = 1:size(P2,2)
           d_ = distPolyPoint(P1,P2(:,i));
           d = max(d,d_);
        end
        
    elseif isa(P2,'mptPolytope') || isa(P2,'interval') || ...
           isa(P2,'zonotope') ||  isa(P2,'conZonotope') || ...
           isa(P2,'zonoBundle')
       
       % convert set to polytope
       P2 = mptPolytope(P2);
       
       % compute distance d(P1,P2) = sup_{x \in P1} inf_{y \in P2) d(x,y)
       V = vertices(P1);
       
       for i = 1:size(V,2)
           d_ = distPolyPoint(P2,V(:,i));
           d = max(d,d_);
       end
       
       % compute distance d(P2,P1) = sup_{x \in P2} inf_{y \in P1) d(x,y)
       V = vertices(P2);
       
       for i = 1:size(V,2)
           d_ = distPolyPoint(P1,V(:,i));
           d = max(d,d_);
       end
       
    else
        % throw error for given arguments
        error(noops(P1,P2));
    end
end


% Auxiliary Functions -----------------------------------------------------

function d = distPolyPoint(poly,x)
% compute Hausdorff distance between a polytope and a single point
% according to Equation (7) in [1]

    H = 2*eye(length(x));
    f = -2*x;
    
    A = poly.P.A; 
    b = poly.P.b;
    
    options = optimoptions('quadprog','Display','off');
    
    [~,d] = quadprog(H,f,A,b,[],[],[],[],[],options);

    d = sqrt(d + x'*x);
end

%------------- END OF CODE --------------