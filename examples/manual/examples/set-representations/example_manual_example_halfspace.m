function example_manual_example_halfspace()
% example_manual_example_halfspace - example from the manual demonstrating 
% the halfspace example from the manual
%
% Syntax:
%   example_manual_example_halfspace()
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

% Authors:        Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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
poly = polytope([1 0;-1 0;0 1;0 -1;1 1],[3;1;2;2;2]);

poly_ = H & poly;

plot(poly_,[1,2],'FaceColor',[0 .7 0]);
plot(poly,[1,2],'b');

% plot --------------------------------------------------------------------

% enlargeAxis(1.2)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
