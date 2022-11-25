function R = guardIntersect_polytope(obj,R,guard,options)
% guardIntersect_polytope - enclosure of guard intersections based on
%                           using a combination of zonotopes and polytopes
%                           as described in [1]
%
% Syntax:  
%    R = guardIntersect_polytope(obj,R,options)
%
% Inputs:
%    obj - object of class location
%    R - list of intersections between the reachable set and the guard
%    guard - guard set
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - enclosure of the guard intersection
%
% References: 
%   [1] M. Althoff et al. "Computing Reachable Sets of Hybrid Systems Using 
%       a Combination of Zonotopes and Polytopes", 2009
%   [2] M. Althoff et al. "Zonotope bundles for the efficient computation 
%       of reachable sets", 2011
% 
% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      19-December-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % enclose all relevant reachable sets by polytopes
    for i = 1:length(R)
        R{i} = conv2polytope(R{i});
    end

    % intersect the reachable sets with the guard set
    for i = 1:length(R)
       R{i} = R{i} & guard; 
    end

    % compute vertices
    V = [];
    
    for i = 1:length(R)
        if ~isempty(R{i})
            vert = vertices(R{i});
            V = [V,vert];
        end
    end
    
    % enclose vertices with the methods described in Section V.A in [2]
    m = length(options.enclose);
    Z = cell(m,1);
    
    for i = 1:m
        
        switch options.enclose{i}
            case 'box'
                % Box method as described in Section V.A.a) in [2]
                Z{i} = zonotope(interval.enclosePoints(V));
                
            case 'pca'
                % PCA mehtod as described in Section V.A.b) in [2]
                Z{i} = zonotope.enclosePoints(V,'stursberg');
                
            case 'flow'
                % Flow method as described in Section V.A.d) in [2]
                Z{i} = flowEnclosure(obj.contDynamics,V,options); 
                
            otherwise
                error('Wrong value for options.enclose!');
                
        end
    end
    
    % construct the enclosing zonotope bundle object
    if length(Z) == 1
       R = Z{1}; 
    else
       R = zonoBundle(Z); 
    end 
end


% Auxiliary Functions -----------------------------------------------------

function Z = flowEnclosure(sys,V,options)
% Flow method for enclosing a set of vertices with a zonotope as described 
% in Section V.A.d) in [1]

    % compute flow at center of vertices
    c = mean(V,2);
    
    options.u = options.uTrans;
    fcnHan = getfcn(sys,options);
    dir = fcnHan(0,c);
    
    dir = dir./norm(dir);
    
    % get basis that is orthogonal to the flow direction
    B = gramSchmidt(dir);
    
    % compute box enclosure in transformed space 
    Z = B*zonotope(interval.enclosePoints(B'*V));

end

function P = conv2polytope(R)
% compute an enclosing polytope for the given set

    % enclose set with zonotope
    if ~isa(R,'zonotope')
       R = zonotope(R); 
    end

    % reduce the given set with different methods
    Zred1 = reduce(R,'pca',1);
    Zred2 = zonotope(interval(R));
    
    % compute intersection of the two reduced sets
    P = mptPolytope(Zred1) & mptPolytope(Zred2);
end

%------------- END OF CODE --------------