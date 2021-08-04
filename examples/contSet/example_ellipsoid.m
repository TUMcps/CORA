function completed = example_ellipsoid()
% example_ellipsoid - example demonstrating set based computation with
%                     ellispoids
%
% Syntax:  
%    completed = example_ellipsoid()
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

E1 = ellipsoid(diag([1/2,2])) % create ellipsoid E1 and display it
A = diag([2,0.5]);

E2 = A*E1 + 0.5; % linear Map + Minkowski addition
E3 = E1 + E2; % Minkowski addition
E4 = E1 & E2; % intersection

disp(['E1 in E2?: ',num2str(in(E2,E1))]);
disp(['E1 in E3?: ',num2str(in(E3,E1))]);

figure; hold on
plot(E1,[1,2],'b'); % plot E1 in blue
plot(E2,[1,2],'g'); % plot E2 in green
plot(E3,[1,2],'r'); % plot E3 in red
plot(E4,[1,2],'k'); % plot E4 in black

E5 = ellipsoid([0.8,-0.6; -0.6,0.8],[1; -4]); % create ellipsoid E5
Zo_box = zonotope(E5); % over-approximate E5 by a parallelotope
Zu_norm = zonotope(E5,10,'i:norm'); % inner-approximate E5 using zonotope norm

figure; hold on
plot(E5); % plot E5
plot(Zo_box,[1,2],'r');  % plot over-approximative zonotope Zo_box
plot(Zu_norm,[1,2],'m'); % plot under-approximative zonotope Zu_norm

% example completed
completed = 1;

%------------- END OF CODE --------------