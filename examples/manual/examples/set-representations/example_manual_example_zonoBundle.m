function example_manual_example_zonoBundle()
% example_manual_example_zonoBundle - example from the manual demonstrating 
% the zonoBundle example from the manual
%
% Syntax:
%   example_manual_example_zonoBundle()
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

Z{1} = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1;
Z{2} = zonotope([-1 1 0; 1 0 1]); % create zonotope Z2;
Zb = zonoBundle(Z); % instantiate zonotope bundle from Z1, Z2
 
vol = volume(Zb) % compute and display volume of zonotope bundle
figure; hold on
plot(Z{1}); % plot Z1
plot(Z{2}); % plot Z2
plot(Zb,[1 2],'FaceColor',[.675 .675 .675]); % plot Zb in gray

% plot --------------------------------------------------------------------

enlargeAxis(1.2)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
