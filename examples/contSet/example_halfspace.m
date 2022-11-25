function completed = example_halfspace()
% example_halfspace - example instantiation of halfspace objects
%
% Syntax:  
%    completed = example_halfspace()
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author:        ---
% Written:       ---
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% construct halfspace object
c = [1;1];
d = 1;

H = halfspace(c,d);

% visualize the halfspace
figure
hold on
xlim([-2,4]);
ylim([-3,3]);

plot(H,[1,2],'r','FaceAlpha',0.5);

% intersect halfspace with polytope
poly = mptPolytope([1 0;-1 0;0 1;0 -1;1 1],[3;1;2;2;2]);

poly_ = H & poly;

plot(poly_,[1,2],'FaceColor',[0 .7 0],'Filled',true,'EdgeColor','none');
plot(poly,[1,2],'b');


% example completed
completed = 1;

%------------- END OF CODE --------------