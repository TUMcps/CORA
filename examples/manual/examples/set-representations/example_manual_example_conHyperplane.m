function example_manual_example_conHyperplane()
% example_manual_example_conHyperplane - example from the manual demonstrating 
% the conHyperplane example from the manual
%
% Syntax:
%   example_manual_example_conHyperplane()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Tobias Ladner
% Written:       27-September-2023
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
figure
hold on
xlim([-2,4]);
ylim([-3,3]);
plot(conHyperplane(c,d),[1,2],'r');             % unconstrained hyperplane
plot(polytope(A,b),[1,2],'g');               % inequality constraints
plot(hyp,[1,2],'b');                            % constrained hyperplane

% plot --------------------------------------------------------------------

% enlargeAxis(1.2)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
