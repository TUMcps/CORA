function completed = example_probZonotope()
% updated: 21-April-2018, MA

Z1=[10 ; 0 ]; % uncertain center
Z2=[0.6 1.2  ; 0.6 -1.2]; % generators with normally distributed factors
pZ=probZonotope(Z1,Z2,2); % probabilistic zonotope

M=[-1 -1;1 -1]*0.2; % mapping matrix
pZencl = enclose(pZ,M); % probabilistic enclosure of pZ and M*pZ

figure('renderer','zbuffer')
hold on
plot(pZ,'dark'); % plot pZ
plot(expm(M)*pZ,'light'); % plot expm(M)*pZ
plot(pZencl,'mesh') % plot enclosure

campos([-3,-51,1]); %set camera position
drawnow; % draw 3D view

%example completed
completed = 1;

%------------- END OF CODE --------------