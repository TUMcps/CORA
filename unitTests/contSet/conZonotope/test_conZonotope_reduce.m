function res = test_conZonotope_reduce
% test_conZonotope_reduce - unit test function for order reduction of
%                           constrained zonotopes
%
% Syntax:  
%    res = test_conZonotope_reduce
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

methods = {'girard','combastel'};

for j = 1:length(methods)
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
    V = get(v,'V');
    
    % Reduce all constraints for which the elimination does not result in
    % an over-approximation
    cZono = reduce(cZono,'redConstr');

    % Reduce the constrained zonotope
    nc = size(cZono.A,1);
    ng = size(cZono.Z,2)-1;
    n = size(cZono.Z,1);
    o = floor((ng-nc)/n)+1;
    
    cRed{1} = reduce(cZono,methods{j},o,nc-1);   % reduce 1 constraint
    cRed{2} = reduce(cZono,methods{j},o,nc-4);   % reduce 2 constraints
    cRed{3} = reduce(cZono,methods{j},o-1,nc);   % reduce 1 generator
    cRed{4} = reduce(cZono,methods{j},o-4,nc);   % reduce 2 generators
    cRed{5} = reduce(cZono,methods{j},o-1,nc-1); % reduce constraint and generator
    
    % Chech if the vertices of the original constrained zonotope are
    % located inside the reduced constrained zonotope
    for i = 1:length(cRed)
       
        % convert to polytope (for easy checks if point is inside set)
        poly = mptPolytope(cRed{i});
        P = get(poly,'P');
        
%         % plot the result
%         hold on
%         plot(cZono,[1,2],'r');
%         plot(cRed{i},[1,2],'b');
        
        % check if all vertices are located inside the set
        temp = P.A*V-P.b*ones(1,size(V,2));
        if any(any(temp > 1e-10))
            error('Test 1 (test case %i) failed!',i); 
        end
    end
end



% TEST 2: 3D --------------------------------------------------------------

methods = {'girard','combastel'};

for j = 1:length(methods)
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
    V = get(v,'V');
    
    % Reduce all constraints for which the elimination does not result in
    % an over-approximation
    cZono = reduce(cZono,'redConstr');

    % Reduce the constrained zonotope
    nc = size(cZono.A,1);
    ng = size(cZono.Z,2)-1;
    n = size(cZono.Z,1);
    o = floor((ng-nc)/n)+1;
    
    cRed{1} = reduce(cZono,methods{j},o,nc-1);   % reduce 1 constraint
    cRed{2} = reduce(cZono,methods{j},o,nc-4);   % reduce 2 constraints
    cRed{3} = reduce(cZono,methods{j},o-1,nc);   % reduce 1 generator
    cRed{4} = reduce(cZono,methods{j},o-4,nc);   % reduce 2 generators
    cRed{5} = reduce(cZono,methods{j},o-1,nc-1); % reduce constraint and generator
    
    % Chech if the vertices of the original constrained zonotope are
    % located inside the reduced constrained zonotope
    for i = 1:length(cRed)
       
        % convert to polytope (for easy checks if point is inside set)
        poly = mptPolytope(cRed{i});
        P = get(poly,'P');
        
        % check if all vertices are located inside the set
        temp = P.A*V-P.b*ones(1,size(V,2));
        if any(any(temp > 1e-10))
            error('Test 2 (test case %i) failed!',i); 
        end
    end
end


% Test 3: Redundant Constraints (analytical test) -------------------------

load('conZonotope_reduce.mat');

% remove all redundant constraints
cZred = reduce(cZ,'redConstr');

% calculate vertices
temp = vertices(cZred);
V_ = get(temp,'V');
V_ = V_';

% calculate vertices of the equivalent polytope
poly = mptPolytope(cZ);
temp = vertices(poly);
V = get(temp,'V');
V = V';

% compare with the real vertices
for i = 1:size(V_,2)
   if ~ismembertol(V_(i,:),V,'ByRows',true)
       error('Test 3 (remove redundant constraints)') 
   end
end

res = 1;

%------------- END OF CODE --------------