function completed = example_probZonotope()
% example_probZonotope - example instantiation of probZonotope objects
%
% Syntax:
%    completed = example_probZonotope()
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

% Authors:       Matthias Althoff
% Written:       21-April-2018
% Last update:   10-August-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Z1 = [10; 0]; % uncertain center
Z2 = [0.6 1.2; 0.6 -1.2]; % generators with normally distributed factors
pZ = probZonotope(Z1,Z2); % probabilistic zonotope

M = [-1 -1;1 -1]*0.2; % mapping matrix
pZencl = enclose(pZ,M); % probabilistic enclosure of pZ and M*pZ

figure % initialize figure
hold on
camlight headlight
    
plot(pZ,[1 2],'FaceColor',[0.2 0.2 0.2],...
    'EdgeColor','none', 'FaceLighting','phong'); % plot pZ 
    
plot(expm(M)*pZ,[1,2],'FaceColor',[0.5 0.5 0.5],...
    'EdgeColor','none', 'FaceLighting','phong'); % plot expm(M)*pZ

plot(pZencl,[1,2],'k','FaceColor','none') % plot enclosure

campos([-3,-51,1]); % set camera position
drawnow; % draw 3D view

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
