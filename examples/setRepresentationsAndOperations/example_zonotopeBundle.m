function completed = example_zonotopeBundle()
% updated: 21-April-2018, MA

Z{1} = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1;
Z{2} = zonotope([-1 1 0; 1 0 1]); % create zonotope Z2;
Zb = zonotopeBundle(Z); % instantiate zonotope bundle from Z1, Z2
vol = volume(Zb) % compute and display volume of zonotope bundle

figure; hold on
plot(Z{1}); % plot Z1 
plot(Z{2}); % plot Z2 
plotFilled(Zb,[1 2],[.675 .675 .675],'EdgeColor','none'); % plot Zb in gray

%example completed
completed = 1;

%------------- END OF CODE --------------