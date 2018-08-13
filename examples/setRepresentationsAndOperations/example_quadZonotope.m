function completed = example_quadZonotope()
% updated: 21-April-2018, MA

c = [0;0]; % starting point
E1 = diag([-1,0.5]); % generators of factors with identical indices
E2 = [1 1; 0.5 0.3]; % generators of factors with identical indices
F = [-0.5; 1]; % generators of factors with different indices
G = [0.3; 0.3]; % independent generators

qZ = quadZonotope(c,E1,E2,F,G); % instantiate quadratic zonotope
Z = zonotope(qZ) % over-approximate by a zonotope

figure; hold on
plot(Z); % plot Z 
plotFilled(qZ,[1 2],7,[],[.6 .6 .6],'EdgeColor','none'); % plot qZ 

%example completed
completed = 1;

%------------- END OF CODE --------------