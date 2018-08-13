function completed = example_vertices()
% updated: 21-April-2018, MA

Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1
V1 = vertices(Z1); % compute vertices of Z1
A = [0.5 1; 1 0.5]; % numerical matrix A

V2{1} = A*V1; % linear map of vertices
V2{2} = V2{1} + [1; 0]; % translation of vertices
V3 = collect(V2{1},V2); % collect vertices of cell array V2
Zencl = zonotope(V3); % obtain parallelotope containing all vertices

figure
hold on
plot(V2{1},'k+'); % plot V2{1}
plot(V2{2},'ko'); % plot V2{2}
plot(Zencl); % plot Zencl 

%example completed
completed = 1;

%------------- END OF CODE --------------