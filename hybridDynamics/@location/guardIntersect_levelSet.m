function R = guardIntersect_levelSet(loc,R,guard)
% guardIntersect_levelSet - implementation of the guard intersection 
%    enclosure with polynomial zonotopes as described in [1]
%
% Syntax:
%    R = guardIntersect_levelSet(loc,R,guard,options)
%
% Inputs:
%    loc - location object
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

% Authors:       Niklas Kochdumper
% Written:       07-January-2020 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute enclosing interval for all reachable sets 
I = [];

for i = 1:length(R)
   
    % interval enclosure (see Step 2 in Sec. 3.2 in [1])
    I_ = interval(R{i});
    % domain tightening only possible for guard sets with comparison
    % operator '=='; for instant transitions with conditions, level sets
    % have comparison operator '<=' and thus the domain cannot be
    % tightened; however, instant transitions also only have one set in R
    try
        I_ = tightenDomain(guard,I_);
    end
    
    % compute union of all intervals (see Step 3 in Sec. 3.2 in [1])
    I = I | I_;
end

% compute the intersection of the reachable set with the guard set
% (see Step 4 in Sec. 3.2 in [1])
pZ = polyZonotope(I); 
R = and_(guard,pZ,'approx');

% ------------------------------ END OF CODE ------------------------------
