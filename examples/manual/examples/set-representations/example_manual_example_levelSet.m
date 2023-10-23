function example_manual_example_levelSet()
% example_manual_example_levelSet - example from the manual demonstrating 
% the levelSet example from the manual
%
% Syntax:
%   example_manual_example_levelSet()
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

% construct level sets
syms x y
eq = sin(x) + y;

ls1 = levelSet(eq,[x;y],'==');
ls2 = levelSet(eq,[x;y],'<=');

% visualize the level sets
figure;
subplot(1,2,1)
xlim([-1.5,1.5]);
ylim([-1,1]);
plot(ls1,[1,2],'b');

subplot(1,2,2)
xlim([-1.5,1.5]);
ylim([-1,1]);
plot(ls2,[1,2],'Color',[0.9451 0.5529 0.5686]);

% plot --------------------------------------------------------------------

subplot(1,2,1)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

subplot(1,2,2)
xlabel('$x_{(1)}$','Interpreter','latex')
% ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
