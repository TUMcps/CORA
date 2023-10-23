function completed = example_conHyperplane()
% example_conHyperplane - example instantiation of conHyperplane objects
%
% Syntax:
%    completed = example_conHyperplane()
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

% construct constrained hyperplane
c = [1;1];
d = 1;
A = [1 0;-1 0;0 1;0 -1;1 1];
b = [3;1;2;2;2];

hyp = conHyperplane(c,d,A,b);

% visualize the constrained hyperplane
figure; hold on
xlim([-2,4]); ylim([-3,3]);

% unconstrained hyperplane
plot(conHyperplane(c,d));
% inequality constraints
plot(polytope(A,b),[1,2],'FaceColor',colorblind('gray'));
% constrained hyperplane
plot(hyp,[1,2],'Color',colorblind('r'));

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
