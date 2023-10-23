function R = guardIntersect_conZonotope(loc,R,guard,options)
% guardIntersect_conZonotope - constrained zonotope based enclosure of 
%    guard intersections
%
% Syntax:
%    R = guardIntersect_conZonotope(loc,R,guard,options)
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

% convert all relevant reachable sets to constrained zonotopes
for i=1:length(R)
    
    if isa(R{i},'polyZonotope')
        R{i} = zonotope(R{i}); 
    end

    temp = reduce(R{i},options.reductionTechnique,options.guardOrder);
    R{i} = conZonotope(temp);  
end

% intersect the reachable sets with the guard set    
for i = 1:length(R)
    R{i} = and_(R{i},guard,'exact'); 
end

R = R(~cellfun('isempty',R));

% calculate orthogonal basis with the methods in Sec. V.A in [1]
B = calcBasis(loc,R,guard,options);

% loop over all calculated basis 
Z = cell(length(B),1);

for i = 1:length(B)
    
    I = [];
    
    % loop over all reachable sets
    for j = 1:length(R)
       
        % interval enclosure in the transformed space
        intTemp = interval(B{i}'*R{j});
        
        % unite all intervals
        I = I | intTemp;
    end
    
    % backtransformation to the original space
    Z{i} = B{i} * zonotope(I);
end

% construct the enclosing zonotope bundle object
if length(Z) == 1
    R = Z{1}; 
else
    R = zonoBundle(Z); 
end 

% ------------------------------ END OF CODE ------------------------------
