function res = test_conZonotope_polytope
% test_conZonotope_polytope - unit test function for conversion between
%                             constrained zonotopes and polytopes
%
% Syntax:  
%    res = test_conZonotope_polytope
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
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

res = 0;

% TEST 1: 2D --------------------------------------------------------------

for j = 1:5
    % Generate random polytope vertices
    points = rand(2,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    % Construct a mptPolytope object from the vertices
    poly = mptPolytope(V');

    % Convert to constrained zonotope object
    cZono = conZonotope(poly);

    % Calculate vertices
    v = vertices(cZono);
    V1 = get(v,'V');

    % Convert back to a mptPolytope
    poly = mptPolytope(cZono);

    % Calculate vertices
    v = vertices(poly);
    V2 = get(v,'V');

    % plot the result
%     plotFilled(cZono,[1,2],'b','EdgeColor','none');
%     hold on
%     plot(poly,[1,2],'r');
%     plot(V(1,:),V(2,:),'.k','MarkerSize',12);


    % Check for correctness
    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V1',1e-10,'ByRows',true)
          error('Test 1 failed! Wrong conZonotope vertices!'); 
       end
    end

    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V2',1e-10,'ByRows',true)
          error('Test 1 failed! Wrong mptPolytope vertices!'); 
       end
    end
end


% TEST 2: 3D --------------------------------------------------------------

for j = 1:5
    % Generate random polytope vertices
    points = rand(3,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    % Construct a mptPolytope object from the vertices
    poly = mptPolytope(V');

    % Convert to constrained zonotope object
    cZono = conZonotope(poly);

    % Calculate vertices
    v = vertices(cZono);
    V1 = get(v,'V');

    % Convert back to a mptPolytope
    poly = mptPolytope(cZono);

    % Calculate vertices
    v = vertices(poly);
    V2 = get(v,'V');

    % Check for correctness
    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V1',1e-10,'ByRows',true)
          error('Test 2 failed! Wrong conZonotope vertices!'); 
       end
    end

    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V2',1e-10,'ByRows',true)
          error('Test 2 failed! Wrong mptPolytope vertices!'); 
       end
    end
end




res = 1;


%------------- END OF CODE --------------