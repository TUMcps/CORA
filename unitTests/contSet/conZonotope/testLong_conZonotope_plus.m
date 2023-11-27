function res = testLong_conZonotope_plus
% testLong_conZonotope_plus - unit test function for the
%    calculation of the Minkowski sum of a constrained zonotope object with
%    another set
%
% Syntax:
%    res = testLong_conZonotope_plus
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

% Authors:       Niklas Kochdumper
% Written:       28-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% TEST 1: Random Test (zonotope 2D) ---------------------------------------

% Generate random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

P = polytope(V);
cZ = conZonotope(P);

% generate random zonotope object
Z = zonotope(rand(2,5)-0.5*ones(2,5));

% Minkowski sum of constrained zonotope and zonotope object
cZ_ = cZ + Z;

% calculate points that have to be located inside the resulting conZonotope
VZ = vertices(Z);

N = size(V,2) * size(VZ,2);
points = zeros(2,N);

counter = 1;

for i = 1:size(V,2)
    for j = 1:size(VZ,2)
        points(:,counter) = V(:,i) + VZ(:,j);
        counter = counter + 1;
    end
end

% convert the resulting conZonotope to a polytope (to easily check if
% a point is located inside the conZonotope)
P = polytope(cZ_);

% extract inequality constraints
A = P.A;
b = P.b;

% % visualize result
% plot(points(1,:),points(2,:),'.k');
% hold on
% plot(cZ_,[1,2],'r');

% check if all points are located inside the resulting conZonotope
for i = 1:size(points,2)
   if ~all(A*points(:,i) < b | withinTol(A*points(:,i),b))
       throw(CORAerror('CORA:testFailed'));
   end
end


% TEST 2: Random Test (interval 2D) ---------------------------------------

% Generate random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

P = polytope(V);
cZ = conZonotope(P);

% generate random interval object
temp1 = rand(2,1)-0.5*ones(2,1);
temp2 = rand(2,1)-0.5*ones(2,1);

inter = interval(min(temp1,temp2),max(temp1,temp2));

% Minkowski sum of constrained zonotope and zonotope object
cZ_ = cZ + inter;

% calculate points that have to be located inside the resulting conZonotope
Vinter = vertices(inter);

N = size(V,2) * size(Vinter,2);
points = zeros(2,N);

counter = 1;

for i = 1:size(V,2)
    for j = 1:size(Vinter,2)
        points(:,counter) = V(:,i) + Vinter(:,j);
        counter = counter + 1;
    end
end

% convert the resulting conZonotope to a polytope (to easily check if
% a point is located inside the conZonotope)
P = polytope(cZ_);

% extract inequality constraints
A = P.A;
b = P.b;

% % visualize result
% plot(points(1,:),points(2,:),'.k');
% hold on
% plot(cZ_,[1,2],'r');

% check if all points are located inside the resulting conZonotope
for i = 1:size(points,2)
   if ~all(A*points(:,i) < b | withinTol(A*points(:,i),b))
       throw(CORAerror('CORA:testFailed'));
   end
end


% TEST 3: Random Test (conZonotope 2D) ------------------------------------

% Generate random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

P = polytope(V);
cZ = conZonotope(P);

% generate a second random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V2 = points(:,ind);

P = polytope(V2);
cZ2 = conZonotope(P);

% Minkowski sum of constrained zonotope and zonotope object
cZ_ = cZ + cZ2;

% calculate points that have to be located inside the resulting conZonotope
N = size(V,2) * size(V2,2);
points = zeros(2,N);

counter = 1;

for i = 1:size(V,2)
    for j = 1:size(V2,2)
        points(:,counter) = V(:,i) + V2(:,j);
        counter = counter + 1;
    end
end

% convert the resulting conZonotope to a polytope (to easily check if
% a point is located inside the conZonotope)
P = polytope(cZ_);

% extract inequality constraints
A = P.A;
b = P.b;

% % visualize result
% plot(points(1,:),points(2,:),'.k');
% hold on
% plot(cZ_,[1,2],'r');

% check if all points are located inside the resulting conZonotope
for i = 1:size(points,2)
   if ~all(A*points(:,i) < b | withinTol(A*points(:,i),b,1e-6))
       throw(CORAerror('CORA:testFailed'));
   end
end


% TEST 4: Random Test (vector 2D) -----------------------------------------

% Generate random conZonotope object
points = rand(2,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

P = polytope(V);
cZ = conZonotope(P);

% generate a random vector
vec = rand(2,1) - 0.5*ones(2,1);

% Minkowski sum of constrained zonotope and zonotope object
cZ_ = cZ + vec;

% calculate points that have to be located inside the resulting conZonotope
for i = 1:size(points,2)
   points(:,i) = points(:,i) + vec; 
end

% convert the resulting conZonotope to a polytope (to easily check if
% a point is located inside the conZonotope)
P = polytope(cZ_);

% extract inequality constraints
A = P.A;
b = P.b;

% % visualize result
% plot(points(1,:),points(2,:),'.k');
% hold on
% plot(cZ_,[1,2],'r');

% check if all points are located inside the resulting conZonotope
Tol = 1e-12;
for i = 1:size(points,2)
   if any(A*points(:,i) - b > Tol)
       throw(CORAerror('CORA:testFailed'));
   end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
