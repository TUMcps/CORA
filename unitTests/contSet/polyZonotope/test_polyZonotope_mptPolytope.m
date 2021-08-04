function res = test_polyZonotope_mptPolytope
% test_polyZonotope_mptPolytope - unit test function for the conversion of
%    a mptPolytope to a polynomial zonotope
%
% Syntax:  
%    res = test_polyZonotope_mptPolytope
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

% Author:       Niklas Kochdumper
% Written:      30-October-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 0;

%% ANALYTICAL TEST

% construct a polytope by computation of the convex hull between a zonotope
% and one point
pZ1 = polyZonotope([0;0],[1 0 -1;1 1 1],[],eye(3));
pZ2 = polyZonotope([1.5;2.5],[],[],[]);

pZ = enclose(pZ1,pZ2);

% convert the polynomial zonotope to a mptPolytope
poly = mptPolytope(pZ);

% % visualize the result
% hold on
% plot(pZ,[1,2],'r','Splits',12,'EdgeColor','none');
% plot(poly,[1,2],'b');

% calculate the vertices
V = vertices(poly);

% define the ground truth
V_ = [-2 -2 0 1.5 2 2 0;-1 1 3 2.5 1 -1 -3];

% check for correctness
for i = 1:size(V,2)
   if ~ismembertol(V(:,i)',V_',1e-12,'ByRows',true) 
      error('test_polyZonotope_mptPolytope: analytical test failed!');
   end
end





res = 1;

%------------- END OF CODE --------------