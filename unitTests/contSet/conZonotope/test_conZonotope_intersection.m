function res = test_conZonotope_intersection
% test_conZonotope_intersection - unit test function for intersection of a
%                                 constrained zonotope with other sets
%
% Syntax:  
%    res = test_conZonotope_intersection
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
v = vertices(intZono);
V = get(v,'V');

% define ground truth
V_ = [1 3 1 3;11/8 -1/8 -11/12 -7/12];

% % plot the result
% plot(cZono1,[1,2],'r');
% hold on
% plot(cZono2,[1,2],'b');
% plotFilled(intZono,[1,2],'g','EdgeColor','none');
% plot(V_(1,:),V_(2,:),'.k','MarkerSize',12);

% check correctness
for i = 1:size(V_,2)
   if ~ismembertol(V_(:,i)',V','ByRows',true)
      error('Test 1 (conZonotope analytical) failed!'); 
   end
end

res = 1;



% TEST 2: conZonotope (random) --------------------------------------------

% loop over different dimensions
for j = 2:3
    
    % Generate random polytope vertices 1
    points = rand(j,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);
    poly1 = mptPolytope(V');

    % Generate random polytope vertices 1
    points = rand(j,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);
    poly2 = mptPolytope(V');

    % calculate constrained zonontope intersection
    cZono1 = conZonotope(poly1);
    cZono2 = conZonotope(poly2);
    zonoInt = cZono1 & cZono2;
    v = vertices(zonoInt);
    V = get(v,'V');

    % calculate vertices from polytope interesection
    polyInt = poly1 & poly2;
    v = vertices(polyInt);
    V_ = get(v,'V');

    % plot the result
%     if j == 2
%         plot(cZono1,[1,2],'r');
%         hold on
%         plot(cZono2,[1,2],'b');
%         plotFilled(zonoInt,[1,2],'g','EdgeColor','none');
%         plot(V(1,:),V(2,:),'.k','MarkerSize',12);
%     end

    % check correctness
    if ~isempty(V_)
        for i = 1:size(V_,2)
           if ~ismembertol(V_(:,i)',V',1e-10,'ByRows',true)
              error('Test 2 (conZonotope random) failed!'); 
           end
        end
    end
end



% TEST 3: halfspace (analytical) ------------------------------------------

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
v = vertices(intZono);
V = get(v,'V');

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
      error('Test 3 (halfspace analytical) failed!'); 
   end
end




% TEST 4: constrainedHyperplane (analytical) ------------------------------

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
ch = constrainedHyperplane(hp,Ch,dh);

% calculate intersection
intZono = cZono & ch;
v = vertices(intZono);
V = get(v,'V');

% define ground truth
V_ = [2 2.5;0.5 0.75];

% % plot the result
% x = -4:0.1:4;
% y = (d-C(1)*x)./C(2);
% poly = mptPolytope([Ch;0 1],[dh;4]);
% plotFilled(poly,[1,2],'m','EdgeColor','none','FaceAlpha',0.5);
% hold on
% plot(x,y,'g');
% plot(cZono,[1,2],'r');
% plot(intZono,[1,2],'b','LineWidth',2);

% check correctness
for i = 1:size(V_,2)
   if ~ismembertol(V_(:,i)',V','ByRows',true)
      error('Test 4 (constrainedHyperplane analytical) failed!'); 
   end
end




res = 1;


%------------- END OF CODE --------------