function res = test_conZonotope_plus
% test_conZonotope_plus - unit test function for the caclulation of the
%                         minkowski sum of a constrained zonotope object
%                         with another set
%
% Syntax:  
%    res = test_conZonotope_plus
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
% Written:      28-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 0;


% TEST 1: Random Test (zonotope 2D) ---------------------------------------

% Generate random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

poly = mptPolytope(V');
cZono = conZonotope(poly);

% generate random zonotope object
zono = zonotope(rand(2,5)-0.5*ones(2,5));

% Minkowski sum of constrained zonotope and zonotope object
cZonoRes = cZono + zono;

% calculate points that have to be located inside the resuling conZonotope
temp = vertices(zono);
Vzono = get(temp,'V');

N = size(V,2) * size(Vzono,2);
points = zeros(2,N);

counter = 1;

for i = 1:size(V,2)
    for j = 1:size(Vzono,2)
        points(:,counter) = V(:,i) + Vzono(:,j);
        counter = counter + 1;
    end
end

% convert the resulting conZonotope to a mptPolytope (to easily check if
% a point is located inside the conZonotope)
poly = mptPolytope(cZonoRes);

% extract inequality constraints
temp = get(poly,'P');
A = temp.A;
b = temp.b;

% % visualize result
% plot(points(1,:),points(2,:),'.k');
% hold on
% plot(cZonoRes,[1,2],'r');

% check if all points are located inside the resulting conZonotope
Tol = 1e-14;

for i = 1:size(points,2)
   if any(A*points(:,i) - b > Tol)
       error('Test 1 failed!');
   end
end




% TEST 2: Random Test (interval 2D) ---------------------------------------

% Generate random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

poly = mptPolytope(V');
cZono = conZonotope(poly);

% generate random interval object
temp1 = rand(2,1)-0.5*ones(2,1);
temp2 = rand(2,1)-0.5*ones(2,1);

inter = interval(min(temp1,temp2),max(temp1,temp2));

% Minkowski sum of constrained zonotope and zonotope object
cZonoRes = cZono + inter;

% calculate points that have to be located inside the resuling conZonotope
temp = vertices(inter);
Vinter = get(temp,'V');

N = size(V,2) * size(Vinter,2);
points = zeros(2,N);

counter = 1;

for i = 1:size(V,2)
    for j = 1:size(Vinter,2)
        points(:,counter) = V(:,i) + Vinter(:,j);
        counter = counter + 1;
    end
end

% convert the resulting conZonotope to a mptPolytope (to easily check if
% a point is located inside the conZonotope)
poly = mptPolytope(cZonoRes);

% extract inequality constraints
temp = get(poly,'P');
A = temp.A;
b = temp.b;

% % visualize result
% plot(points(1,:),points(2,:),'.k');
% hold on
% plot(cZonoRes,[1,2],'r');

% check if all points are located inside the resulting conZonotope
Tol = 1e-14;

for i = 1:size(points,2)
   if any(A*points(:,i) - b > Tol)
       error('Test 2 failed!');
   end
end





% TEST 3: Random Test (conZonotope 2D) ------------------------------------

% Generate random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

poly = mptPolytope(V');
cZono = conZonotope(poly);

% generate a second random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V2 = points(:,ind);

poly = mptPolytope(V2');
cZono2 = conZonotope(poly);

% Minkowski sum of constrained zonotope and zonotope object
cZonoRes = cZono + cZono2;

% calculate points that have to be located inside the resuling conZonotope
N = size(V,2) * size(V2,2);
points = zeros(2,N);

counter = 1;

for i = 1:size(V,2)
    for j = 1:size(V2,2)
        points(:,counter) = V(:,i) + V2(:,j);
        counter = counter + 1;
    end
end

% convert the resulting conZonotope to a mptPolytope (to easily check if
% a point is located inside the conZonotope)
poly = mptPolytope(cZonoRes);

% extract inequality constraints
temp = get(poly,'P');
A = temp.A;
b = temp.b;

% % visualize result
% plot(points(1,:),points(2,:),'.k');
% hold on
% plot(cZonoRes,[1,2],'r');

% check if all points are located inside the resulting conZonotope
Tol = 1e-14;

for i = 1:size(points,2)
   if any(A*points(:,i) - b > Tol)
       error('Test 3 failed!');
   end
end




% TEST 4: Random Test (vector 2D) -----------------------------------------

% Generate random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

poly = mptPolytope(V');
cZono = conZonotope(poly);

% generate a random vector
vec = rand(2,1) - 0.5*ones(2,1);

% Minkowski sum of constrained zonotope and zonotope object
cZonoRes = cZono + vec;

% calculate points that have to be located inside the resuling conZonotope
for i = 1:size(points,2)
   points(:,i) = points(:,i) + vec; 
end

% convert the resulting conZonotope to a mptPolytope (to easily check if
% a point is located inside the conZonotope)
poly = mptPolytope(cZonoRes);

% extract inequality constraints
temp = get(poly,'P');
A = temp.A;
b = temp.b;

% % visualize result
% plot(points(1,:),points(2,:),'.k');
% hold on
% plot(cZonoRes,[1,2],'r');

% check if all points are located inside the resulting conZonotope
Tol = 1e-14;

for i = 1:size(points,2)
   if any(A*points(:,i) - b > Tol)
       error('Test 4 failed!');
   end
end

res = 1;