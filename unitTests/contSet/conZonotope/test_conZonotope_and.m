function res = test_conZonotope_and
% test_conZonotope_and - unit test function for intersection of a
%                        constrained zonotope with other sets
%
% Syntax:  
%    res = test_conZonotope_and
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

res = false;

% TEST 1: conZonotope (analytical) ----------------------------------------

% constrained zonotope 1
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZono1 = conZonotope(Z,A,b);

% constrained zonotope 2
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZono2 = conZonotope(Z,A,b);

% calculate intersection
intZono = cZono1 & cZono2;
V = vertices(intZono);

% define ground truth
V_ = [1 3 1 3;11/8 -1/8 -11/12 -7/12];

% % plot the result
% plot(cZono1,[1,2],'r');
% hold on
% plot(cZono2,[1,2],'b');
% plot(intZono,[1,2],'g','Filled',true,'EdgeColor','none');
% plot(V_(1,:),V_(2,:),'.k','MarkerSize',12);

% check correctness
for i = 1:size(V_,2)
   if ~ismembertol(V_(:,i)',V','ByRows',true)
      error('Test 1 (conZonotope analytical) failed!'); 
   end
end

% TEST 2: halfspace (analytical) ------------------------------------------

% constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZono = conZonotope(Z,A,b);

% hyperplane
C = [1 -2];
d = 1;
hp = halfspace(C,d);

% calculate intersection
intZono = cZono & hp;
V = vertices(intZono);

% define ground truth
V_ = [1 3;0 1];

% % plot the result
% x = -4:0.1:4;
% y = (d-C(1)*x)./C(2);
% hold on
% plot(x,y,'g');
% plot(cZono,[1,2],'r');
% plot(intZono,[1,2],'b');

% check correctness
for i = 1:size(V_,2)
   if ~ismembertol(V_(:,i)',V','ByRows',true)
      error('Test 2 (halfspace analytical) failed!'); 
   end
end




% TEST 3: conHyperplane (analytical) --------------------------------------

% constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZono = conZonotope(Z,A,b);

% constrained hyperplane
C = [1 -2];
d = 1;
hp = halfspace(C,d);
Ch = [-2 -0.5;1 0];
dh = [-4.25;2.5];
ch = conHyperplane(hp,Ch,dh);

% calculate intersection
intZono = cZono & ch;
V = vertices(intZono);

% define ground truth
V_ = [2 2.5;0.5 0.75];

% % plot the result
% x = -4:0.1:4;
% y = (d-C(1)*x)./C(2);
% poly = mptPolytope([Ch;0 1],[dh;4]);
% plot(poly,[1,2],'m','Filled',true,'EdgeColor','none','FaceAlpha',0.5);
% hold on
% plot(x,y,'g');
% plot(cZono,[1,2],'r');
% plot(intZono,[1,2],'b','LineWidth',2);

% check correctness
for i = 1:size(V_,2)
   if ~ismembertol(V_(:,i)',V','ByRows',true)
      error('Test 3 (conHyperplane analytical) failed!'); 
   end
end


res = true;


%------------- END OF CODE --------------