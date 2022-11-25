function R = guardIntersect_nondetGuard(obj,R,guard,options)
% guardIntersect_nondetGuard - enclosure of guard intersections for 
%                              non-deterministic guard sets with large 
%                              uncertainty
%
% Syntax:  
%    R = guardIntersect_nondetGuard(obj,R,guard,options)
%
% Inputs:
%    obj - object of class location
%    R - list of intersections between the reachable set and the guard
%    guard - guard set (class: constrained hyperplane)
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - set enclosing the guard intersection
%
% References: 
%   [1] M. Althoff et al. "Zonotope bundles for the efficient computation 
%       of reachable sets", 2011
% 
% Author:       Niklas Kochdumper
% Written:      19-December-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % calc. orthogonal basis with the methods described in Sec. V.A in [1]
    B = calcBasis(obj,R,guard,options);
   
    % loop over all basis
    Z = cell(length(B),1);
        
    for i = 1:length(B)
       
        % enclose all reachable set with an interval in transformed space
        int = [];
        
        for j = 1:length(R)
            intnew = interval(B{i}'*R{j});
            if isempty(intnew) && isempty(int)
                int = int | B{i}'*interval(R{j});
            else
                int = int | intnew;
            end
        end

        % backtransformation to the original space
        Ztemp = B{i} * zonotope(int);
        
        % convert the resulting set to a constrained zonotope
        cZ = conZonotope(Ztemp);
        
        % intersect the set with the guard set
        cZ = cZ & guard;
        
        % enclose the intersection with an interval in transformed space
        Z{i} = B{i} * zonotope(interval(B{i}'*cZ));
    end
    
    % construct set enclosing the intersection
    if length(Z) == 1
        R = Z{1}; 
    else
        R = zonoBundle(Z); 
    end
end

%------------- END OF CODE --------------