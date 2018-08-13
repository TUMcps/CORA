function res = test_conZonotope_interval
% test_conZonotope_interval - unit test function for the caclulation of
%                             a bounding box of a constrained zonotope object
%
% Syntax:  
%    res = test_conZonotope_interval
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
% Written:      22-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 0;

% TEST 1: Figure 1 in [1] -------------------------------------------------

% construct zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZono = conZonotope(Z,A,b);

% calculate interval
int = interval(cZono);

% % plot the result
% plot(cZono,[1,2],'r');
% hold on
% plot(int,[1,2],'b');

% compare with ground-truth
int_ = interval([-2.5;-1.5],[3.5;2.5]);

for i = 1:length(int)
   if infimum(int_(i)) ~= infimum(int(i)) || supremum(int_(i)) ~= supremum(int(i))
      error('Test 1 failed!'); 
   end
end



% TEST 2: Figure 2 in [1] -------------------------------------------------

% construct zonotope
Z = [0 1 0 1;0 1 2 -1];
A = [-2 1 -1];
b = 2;
cZono = conZonotope(Z,A,b);

% calculate interval
int = interval(cZono);

% % plot the result
% plot(cZono,[1,2],'r');
% hold on
% plot(int,[1,2],'b');

% compare with ground-truth
int_ = interval([-2;-2],[0;3]);

for i = 1:length(int)
   if infimum(int_(i)) ~= infimum(int(i)) || supremum(int_(i)) ~= supremum(int(i))
      error('Test 2 failed!'); 
   end
end



% TEST 3: Random Test 3D --------------------------------------------------

% Generate random conZonotope object
points = rand(3,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

poly = mptPolytope(V');
cZono = conZonotope(poly);

% calculate interval
int = interval(cZono);

% compare with ground-truth for the vertices
v = vertices(cZono);
V = get(v,'V');
int_ = interval(min(V,[],2),max(V,[],2));

for i = 1:length(int)
   if abs(infimum(int_(i)) - infimum(int(i))) > 1e-10 || abs(supremum(int_(i)) - supremum(int(i))) > 1e-10
      error('Test 3 failed!'); 
   end
end


res = 1;