function example_manual_example_zonotope()
% example_manual_example_zonotope - example from the manual demonstrating 
% the zonotope example from the manual
%
% Syntax:
%   example_manual_example_zonotope()
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

Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1
Z2 = zonotope([-1 1 0; 1 0 1]); % create zonotope Z2
A = [0.5 1; 1 0.5]; % numerical matrix A

Z3 = Z1 + Z2; % Minkowski addition
Z4 = A*Z3; % linear map
 
f1 = figure; hold on
plot(Z1,[1 2],'b'); % plot Z1 in blue
plot(Z2,[1 2],'g'); % plot Z2 in green
plot(Z3,[1 2],'r'); % plot Z3 in red
plot(Z4,[1 2],'k'); % plot Z4 in black

P = polytope(Z4) % convert to and display halfspace representation
I = interval(Z4) % convert to and display interval

f2 = figure; hold on
plot(Z4); % plot Z4
plot(I,[1 2],'g'); % plot intervalhull in green

% plot --------------------------------------------------------------------

figure(f1);
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')
enlargeAxis(1.2)

figure(f2);
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
