function res = testLong_conZonotope_and
% testLong_conZonotope_and - unit test function for intersection
%    of a constrained zonotope with other sets
%
% Syntax:  
%    res = testLong_conZonotope_and
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% TEST 1: conZonotope (random) --------------------------------------------

% loop over different dimensions
for j = 2:3
    
    % Generate random polytope vertices 1
    points = rand(j,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);
    P1 = mptPolytope(V');

    % Generate random polytope vertices 1
    points = rand(j,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);
    P2 = mptPolytope(V');

    % calculate constrained zonontope intersection
    cZ1 = conZonotope(P1);
    cZ2 = conZonotope(P2);
    zonoInt = cZ1 & cZ2;
    V = vertices(zonoInt);

    % calculate vertices from polytope interesection
    polyInt = P1 & P2;
    V_ = vertices(polyInt);

    % plot the result
%     if j == 2
%         plot(cZono1,[1,2],'r');
%         hold on
%         plot(cZono2,[1,2],'b');
%         plot(zonoInt,[1,2],'FaceColor','g');
%         plot(V(1,:),V(2,:),'.k','MarkerSize',12);
%     end

    % check correctness
    if ~compareMatrices(V,V_,1e-10)
          res = false;
          break
    end
end

%------------- END OF CODE --------------