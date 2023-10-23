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
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ---
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct halfspace object
c = [1;1];
d = 1;

hs = halfspace(c,d);

% visualize the halfspace
figure; hold on;
xlim([-2,4]);
ylim([-3,3]);

plot(hs,[1,2],'FaceColor',colorblind('r'));

% intersect halfspace with polytope
poly = polytope([1 0;-1 0;0 1;0 -1;1 1],[3;1;2;2;2]);

poly_ = hs & poly;

plot(poly_,[1,2],'FaceColor',colorblind('gray'));
plot(poly);


% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
