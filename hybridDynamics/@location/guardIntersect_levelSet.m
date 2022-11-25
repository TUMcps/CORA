function R = guardIntersect_levelSet(obj,R,guard)
% guardIntersect_zonoGirard - implementation of the guard intersection 
%                             enclosure with polynomial zonotopes as 
%                             described in [1]
%
% Syntax:  
%    R = guardIntersect_levelSet(obj,R,guard,options)
%
% Inputs:
%    obj - object of class location
%    R - list of intersections between the reachable set and the guard
%    guard - guard set (class: levelSet)
%
% Outputs:
%    R - polynomial zonotope enclosing the guard intersection
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: guardIntersect
%
% References: 
%   [1] N. Kochdumper et al. "Reachability Analysis for Hybrid Systems with 
%       Nonlinear Guard Sets", HSCC 2020

% Author: Niklas Kochdumper
% Written: 07-January-2020 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % compute enclosing interval for all reachable sets 
    int = [];
    
    for i = 1:length(R)
       
        % interval enclosure (see Step 2 in Sec. 3.2 in [1])
        int_ = interval(R{i});
        int_ = tightenDomain(guard,int_);
        
        % compute union of all intervals (see Step 3 in Sec. 3.2 in [1])
        int = int | int_;
    end
    
    % compute the intersection of the reachable set with the guard set
    % (see Step 4 in Sec. 3.2 in [1])
    pZ = polyZonotope(int); 
    R = guard & pZ;   
end

%------------- END OF CODE --------------