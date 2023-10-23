function example_manual_example_ellipsoid()
% example_manual_example_ellipsoid - example from the manual demonstrating 
% the ellipsoid example from the manual
%
% Syntax:
%   example_manual_example_ellipsoid()
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

E1 = ellipsoid(diag([1/2,2])) % create ellipsoid E1 and display it
A = diag([2,0.5]);

E2 = A*E1 + 0.5; % linear Map + Minkowski addition
E3 = E1 + E2; % Minkowski addition
E4 = E1 & E2; % intersection
 
disp(['E1 in E2?: ',num2str(E2.contains(E1))]);
disp(['E1 in E3?: ',num2str(E3.contains(E1))]);

f1 = figure; hold on
plot(E1,[1,2],'b'); % plot E1 in blue
plot(E2,[1,2],'g'); % plot E2 in green
plot(E3,[1,2],'r'); % plot E3 in red
plot(E4,[1,2],'k'); % plot E4 in black

E5 = ellipsoid([0.8,-0.6; -0.6,0.8],[1; -4]); % create ellipsoid E5
Zo_box = zonotope(E5); % overapproximate E5 by a parallelotope
Zu_norm = zonotope(E5,10,'outer:norm'); % overapproximate E5 using zonotope norm

f2 = figure; hold on
plot(E5); % plot E5
plot(Zo_box,[1,2],'r'); % plot overapproximative zonotope Zo_box
plot(Zu_norm,[1,2],'m'); % plot overapproximative zonotope Zu_norm

% plot --------------------------------------------------------------------

figure(f1)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

figure(f2)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
