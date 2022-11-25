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


% construct constrained hyperplane
c = [1;1];
d = 1;
A = [1 0;-1 0;0 1;0 -1;1 1];
b = [3;1;2;2;2];

cH = conHyperplane(c,d,A,b);

% visualize the constrained hyperplane
figure
hold on
xlim([-2,4]);
ylim([-3,3]);

plot(conHyperplane(c,d),[1,2],'r');             % unconstrained hyperplane
plot(mptPolytope(A,b),[1,2],'g');               % inequality constraints

plot(cH,[1,2],'b');                             % constrained hyperplane

% example completed
completed = 1;

%------------- END OF CODE --------------