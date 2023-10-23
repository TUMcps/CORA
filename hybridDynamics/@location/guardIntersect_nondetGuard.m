function R = guardIntersect_nondetGuard(loc,R,guard,options)
% guardIntersect_nondetGuard - enclosure of guard intersections for 
%    non-deterministic guard sets with large uncertainty
%
% Syntax:
%    R = guardIntersect_nondetGuard(loc,R,guard,options)
%
% Inputs:
%    loc - location object
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

% Authors:       Niklas Kochdumper
% Written:       19-December-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% calc. orthogonal basis with the methods described in Sec. V.A in [1]
B = calcBasis(loc,R,guard,options);

% loop over all basis
Z = cell(length(B),1);
    
for i = 1:length(B)
   
    % enclose all reachable set with an interval in transformed space
    I = [];
    
    for j = 1:length(R)
        intnew = interval(B{i}'*R{j});
        if representsa_(intnew,'emptySet',eps) && representsa_(I,'emptySet',eps)
            I = I | B{i}'*interval(R{j});
        else
            I = I | intnew;
        end
    end

    % backtransformation to the original space
    Ztemp = B{i} * zonotope(I);
    
    % convert the resulting set to a constrained zonotope
    cZ = conZonotope(Ztemp);
    
    % intersect the set with the guard set
    cZ = and_(cZ,guard,'exact');
    
    % enclose the intersection with an interval in transformed space
    Z{i} = B{i} * zonotope(interval(B{i}'*cZ));
end

% construct set enclosing the intersection
if length(Z) == 1
    R = Z{1}; 
else
    R = zonoBundle(Z); 
end

% ------------------------------ END OF CODE ------------------------------
