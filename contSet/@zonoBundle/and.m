function zB = and(zB1,zB2)
% and - returns the intersection of two zonotope bundles
%
% Syntax: 
%    zB = and(zB1,zB2)
%
% Inputs:
%    zB1 - zonotope bundle
%    zB2 - zonotope bundle or zonotope
%
% Outputs:
%    zB - zonotope bundle after intersection
%
% Example: 
%    % define sets
%    zono1 = zonotope([0 1 2 0;0 1 0 2]);
%    zono2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({zono1,zono2});
%
%    hs = halfspace([1 1],2);
%
%    % intersection
%    res1 = zB & hs;
%
%    % visualization
%    figure
%    hold on
%    xlim([-1,4]);
%    ylim([-4,4]);
%    plot(hs,[1,2],'r','FaceAlpha',0.5);
%    plot(res1,[1,2],'g','Filled',true,'EdgeColor','none');
%    plot(zB,[1,2],'b','LineWidth',3);
%
% Other m-files required: none
% Subfunctions: ---
% MAT-files required: none
%
% See also: zonotope/and

% Author:       Matthias Althoff
% Written:      16-November-2010 
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% determine zonotope bundle object
if ~isa(zB1,'zonoBundle')
   temp = zB1;
   zB1 = zB2;
   zB2 = temp;
end

% different cases for the different types of objects
if isa(zB2,'zonotope')
    
    zB = zB1;
    zB.Z{end+1} = zB2;
    zB.parallelSets = zB.parallelSets + 1;
    
elseif isa(zB2,'zonoBundle')
    
    zB = zB1;
    
    for i = 1:zB2.parallelSets
        zB.Z{end+1} = zB2.Z{i};
    end
    
    zB.parallelSets = zB.parallelSets + zB2.parallelSets;
    
elseif isa(zB2,'interval')
    
    zB = zB1;
    zB.Z{end+1} = zonotope(zB2);
    zB.parallelSets = zB.parallelSets + 1;
    
elseif isa(zB2,'mptPolytope') || isa(zB2,'conZonotope')
    
    zB = zB1 & zonoBundle(zB2);
    
elseif isa(zB2,'halfspace')
    
    % construct basis orthogonal to halfspace normal vector
    B = gramSchmidt(zB2.c);
    
    % compute enclosing interval in transformed space
    zono_ = B' * zB1.Z{1};
    int_ = interval(zono_);
    
    % consider upper bound applied by halfspace constraint c*x <= d
    infi = infimum(int_);
    sup = supremum(int_);
    
    sup(1) = zB2.d/norm(zB2.c);
    
    int_ = interval(infi,sup);
    
    % backtransformation to orginal space
    zono = B * zonotope(int_);
    
    % intersection
    zB = zB1 & zono;
    
    
elseif isa(zB2,'conHyperplane')
    
    % Part 1: intersection with the hyperplane ----------------------------
    
    % construct basis orthogonal to halfspace normal vector
    B = gramSchmidt(zB2.h.c);
    
    % compute enclosing interval in transformed space
    zono_ = B' * zB1.Z{1};
    int_ = interval(zono_);
    
    % consider upper bound applied by halfspace constraint c*x <= d
    infi = infimum(int_);
    sup = supremum(int_);
    
    temp = zB2.h.d/norm(zB2.h.c);
    sup(1) = temp;
    infi(1) = temp;
    
    int_ = interval(infi,sup);
    
    % backtransformation to orginal space
    zono = B * zonotope(int_);
    
    % intersection
    zB = zB1 & zono;
    
   
    % Part 2: intersection with the constraints ---------------------------
    
    % loop over all constraints
    C = zB2.C;
    d = zB2.d;

    for i = 1:size(C,1)

       % construct halfspace
       hs = halfspace(C(i,:)',d(i));

       % check if set is fully contained in halfspace
       if ~in(hs,zB)

           % intersect set with halfspace
           zB = zB & hs; 
       end
    end
    
elseif isa(zB2,'levelSet') || isa(zB2,'conPolyZono')
    
    zB = zB2 & zB1;
    
else
    
    % throw error for given arguments
    error(noops(zB1,zB2));
    
end


%------------- END OF CODE --------------