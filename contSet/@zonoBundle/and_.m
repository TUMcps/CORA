function zB = and_(zB,S,varargin)
% and_ - returns the intersection of a zonotope bundle and another set
%
% Syntax:
%    zB = and_(zB,S)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object
%
% Outputs:
%    zB - zonotope bundle after intersection
%
% Example: 
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%    hs = halfspace([1 1],2);
%
%    res = zB & hs;
%
%    figure; hold on; xlim([-1,4]); ylim([-4,4]);
%    plot(hs,[1,2],'r','FaceAlpha',0.5);
%    plot(res,[1,2],'FaceColor','g');
%    plot(zB,[1,2],'b','LineWidth',3);
%
% Other m-files required: none
% Subfunctions: ---
% MAT-files required: none
%
% See also: contSet/and, zonotope/and_

% Authors:       Matthias Althoff
% Written:       16-November-2010 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: 27-March-2023 (MW, rename and_)

% ------------------------------ BEGIN CODE -------------------------------

% different cases for the different types of objects
if isa(S,'zonotope')
    
    zB.Z{end+1} = S;
    zB.parallelSets = zB.parallelSets + 1;
    
elseif isa(S,'zonoBundle')
    
    % append to list of parallel sets
    for i = 1:S.parallelSets
        zB.Z{end+1} = S.Z{i};
    end
    
    zB.parallelSets = zB.parallelSets + S.parallelSets;
    
elseif isa(S,'interval')
    
    zB.Z{end+1} = zonotope(S);
    zB.parallelSets = zB.parallelSets + 1;
    
elseif isa(S,'polytope') || isa(S,'conZonotope')
    
    zB = and_(zB,zonoBundle(S),'exact');
    
elseif isa(S,'halfspace')
    
    % construct basis orthogonal to halfspace normal vector
    B = gramSchmidt(S.c);
    
    % compute enclosing interval in transformed space
    Z_ = B' * zB.Z{1};
    I_ = interval(Z_);
    
    % consider upper bound applied by halfspace constraint c*x <= d
    lb = infimum(I_);
    ub = supremum(I_);
    
    ub(1) = S.d/norm(S.c);
    
    I_ = interval(lb,ub);
    
    % backtransformation to orginal space
    Z = B * zonotope(I_);
    
    % intersection
    zB = and_(zB,Z,'exact');
    
    
elseif isa(S,'conHyperplane')
    
    % Part 1: intersection with the hyperplane ----------------------------
    
    % construct basis orthogonal to halfspace normal vector
    B = gramSchmidt(S.a');
    
    % compute enclosing interval in transformed space
    Z_ = B' * zB.Z{1};
    I_ = interval(Z_);
    
    % consider upper bound applied by halfspace constraint c*x <= d
    lb = infimum(I_);
    ub = supremum(I_);
    
    temp = S.b/norm(S.a');
    ub(1) = temp;
    lb(1) = temp;
    
    I_ = interval(lb,ub);
    
    % backtransformation to orginal space
    Z = B * zonotope(I_);
    
    % intersection
    zB = and_(zB,Z,'exact');
    
   
    % Part 2: intersection with the constraints ---------------------------
    
    % loop over all constraints
    C = S.C;
    d = S.d;

    for i = 1:size(C,1)

       % construct halfspace
       hs = halfspace(C(i,:)',d(i));

       % check if set is fully contained in halfspace
       if ~contains_(hs,zB)

           % intersect set with halfspace
           zB = and_(zB,hs,'exact');
       end
    end
    
elseif isa(S,'levelSet') || isa(S,'conPolyZono')
    
    zB = and_(S,zB,'exact');
    
else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',zB,S));
    
end

% ------------------------------ END OF CODE ------------------------------
