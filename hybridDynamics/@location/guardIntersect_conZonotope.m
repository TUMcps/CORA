function R = guardIntersect_conZonotope(obj,R,guard,options)
% guardIntersect_conZonotope - constrained zonotope based enclosure of 
%                              guard intersections
%
% Syntax:  
%    R = guardIntersect_conZonotope(obj,R,guard,options)
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

    % convert all relevant reachable sets to constraind zonotopes
    for i = 1:length(R)
            
        if isa(R{i},'polyZonotope')
           R{i} = zonotope(R{i}); 
        end

        temp = reduce(R{i},options.reductionTechnique,options.guardOrder);
        R{i} = conZonotope(temp);  
    end

    % intersect the reachable sets with the guard set    
    for i = 1:length(R)
       R{i} = R{i} & guard; 
    end
    
    R = R(~cellfun('isempty',R));

    % calculate orthogonal basis with the methods in Sec. V.A in [1]
    B = calcBasis(obj,R,guard,options);
    
    % loop over all calculated basis 
    Z = cell(length(B),1);
    
    for i = 1:length(B)
        
        int = [];
        
        % loop over all reachable sets
        for j = 1:length(R)
           
            % interval enclosure in the transformed space
            intTemp = interval(B{i}'*R{j});
            
            % unite all intervals
            int = int | intTemp;
        end
        
        % backtransformation to the original space
        Z{i} = B{i} * zonotope(int);
    end
    
    % construct the enclosing zonotope bundle object
    if length(Z) == 1
       R = Z{1}; 
    else
       R = zonoBundle(Z); 
    end 
end

%------------- END OF CODE --------------