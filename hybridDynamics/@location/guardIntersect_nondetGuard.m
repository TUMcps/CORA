function R = guardIntersect_nondetGuard(loc,R,guard,B)
% guardIntersect_nondetGuard - enclosure of guard intersections for 
%    non-deterministic guard sets with large uncertainty
%
% Syntax:
%    R = guardIntersect_nondetGuard(loc,R,guard,B)
%
% Inputs:
%    loc - location object
%    R - list of intersections between the reachable set and the guard
%    guard - guard set (class: constrained hyperplane)
%    B - basis
%
% Outputs:
%    R - set enclosing the guard intersection

% Authors:       Niklas Kochdumper
% Written:       19-December-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% loop over all basis
Z = cell(length(B),1);
    
for i = 1:length(B)
   
    % enclose all reachable set with an interval in transformed space
    I = interval.empty(size(B{i},1));
    
    for j = 1:length(R)
        I_new = interval(B{i}'*R{j});
        if representsa_(I_new,'emptySet',eps) && representsa_(I,'emptySet',eps)
            I = I | B{i}'*interval(R{j});
        else
            I = I | I_new;
        end
    end

    % backtransformation to the original space
    Z_trans = B{i} * zonotope(I);
    
    % convert the resulting set to a constrained zonotope
    cZ = conZonotope(Z_trans);
    
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
