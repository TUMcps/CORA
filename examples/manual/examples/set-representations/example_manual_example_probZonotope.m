function example_manual_example_probZonotope()
% example_manual_example_probZonotope - example from the manual demonstrating 
% the probZonotope example from the manual
%
% Syntax:
%   example_manual_example_probZonotope()
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

Z1=[10; 0]; % uncertain center
Z2=[0.6 1.2  ; 0.6 -1.2]; % generators with normally distributed factors
pZ=probZonotope(Z1,Z2); % probabilistic zonotope


M=[-1 -1;1 -1]*0.2; % mapping matrix
pZencl = enclose(pZ,M); % probabilistic enclosure of pZ and M*pZ

figure % initialize figure
hold on
camlight headlight

plot(pZ,[1 2],'FaceColor',[0.2 0.2 0.2],...
    'EdgeColor','none', 'FaceLighting','phong'); % plot pZ
plot(expm(M)*pZ,[1,2],'FaceColor',[0.5 0.5 0.5],...
    'EdgeColor','none', 'FaceLighting','phong'); % plot expm(M)*pZ


plot(pZencl,[1,2],'k','FaceAlpha',0) % plot enclosure
campos([-3,-51,1]); % set camera position
drawnow; % draw 3D view

% plot --------------------------------------------------------------------

grid on
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
