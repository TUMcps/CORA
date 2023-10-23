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

E1 = ellipsoid(diag([1/2,2])) % create ellipsoid E1 and display it
A = diag([2,0.5]);

E2 = A*E1 + 0.5; % linear Map + Minkowski addition
E3 = E1 + E2; % Minkowski addition
E4 = E1 & E2; % intersection

disp(['E1 in E2?: ',num2str(contains(E2,E1))]);
disp(['E1 in E3?: ',num2str(contains(E3,E1))]);

figure;

subplot(1,2,1); hold on
% plot original ellipsoid
plot(E1,[1,2],'FaceColor',colorblind('gray'));
% plot linear map and scalar translation
plot(E2);
% plot Minkowksi addition
plot(E3,[1,2],'Color',colorblind('r'));
% plot intersection
plot(E4,[1,2],'k');

E5 = ellipsoid([0.8,-0.6; -0.6,0.8],[1; -4]); % create ellipsoid E5
Zo_box = zonotope(E5); % over-approximate E5 by a parallelotope
Zu_norm = zonotope(E5,10,'inner:norm'); % inner-approximate E5 using zonotope norm

subplot(1,2,2); hold on;
% plot over-approximative zonotope Zo_box
plot(Zo_box,[1,2],'Color',colorblind('r'));
% plot inner-approximative zonotope Zu_norm
plot(Zu_norm,[1,2],'k');
plot(E5);

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
