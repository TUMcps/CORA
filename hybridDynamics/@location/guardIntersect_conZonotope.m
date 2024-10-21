function R = guardIntersect_conZonotope(loc,R,guard,B,options)
% guardIntersect_conZonotope - constrained zonotope based enclosure of 
%    guard intersections
%
% Syntax:
%    R = guardIntersect_conZonotope(loc,R,guard,B,options)
%
% Inputs:
%    loc - location object
%    R - list of intersections between the reachable set and the guard
%    guard - guard set (class: constrained hyperplane)
%    B - basis
%    options - required algorithm parameters: .reductionTechnique, .guardOrder
%
% Outputs:
%    R - set enclosing the guard intersection

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

    R_reduce = reduce(R{i},options.reductionTechnique,options.guardOrder);
    R{i} = conZonotope(R_reduce);  
end

% intersect the reachable sets with the guard set    
for i = 1:length(R)
    R{i} = and_(R{i},guard,'exact'); 
end

R = R(~cellfun('isempty',R));

% loop over all calculated basis 
Z = cell(length(B),1);

for i = 1:length(B)
    
    I = interval.empty(size(B{i},1));
    
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
