function res = testLong_zonotope_reduce
% testLong_zonotope_reduce - unit test function of reduce
%
% Syntax:
%    res = testLong_zonotope_reduce
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       26-July-2016
% Last update:   27-June-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% create zonotope
c = [-1.3; 2.9; 4-2; 3.2];
G = [ ...
0.243, 0.319, 0.022, -0.276, 0.088, -0.019, -0.223, 0.335, -0.517, -0.521, -0.104, 0.659, -0.369, -0.051, -0.746, 0.093, -0.247, 0.119, 0.206, -0.030; ...
-0.436, 0.675, -0.060, -0.209, -0.204, -0.045, -0.008, 0.340, -0.474, 0.252, 0.003, 0.037, 0.817, 0.235, 0.256, 0.019, 0.214, -0.646, 0.739, -0.693; ...
0.259, 0.474, 0.067, -0.760, 0.370, 0.060, -0.042, -0.265, -0.003, -0.153, 0.053, -0.637, -0.214, -0.202, -0.093, 0.013, 0.256, -0.636, -0.191, -0.157; ...
0.589, -0.198, -0.087, 0.371, -0.462, 0.059, 0.162, 0.024, 0.652, 0.757, 0.105, -0.317, 0.260, -0.370, -0.103, 0.104, -0.075, 0.058, 0.054, -0.644];

Zorig = zonotope(c,G);
Zorig_2D = zonotope(c(1:2,:),G(1:2,1:10));


% init result vector
res = [];

%% cluster - order 1
method = 1;
Zred = reduce(Zorig,'cluster',1,[],method);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, 2.617, 3.129, 3.557, 2.983; ...
2.900, -4.435, -0.073 3.883, -2.582; ...
2.000, -4.759, 3.768, -2.589, 1.462; ...
3.200, 0.882, -4.621, -1.789, 2.607];
   
% check result
res(end+1) = compareMatrices(Zmat,true_mat);

% plot
% figure; hold on; box on;
% plot(Zorig); plot(Zred,[1,2],'r'); plot(zonotope(true_mat),[1,2],'g');
% close;

%% cluster - order 3
Zred = reduce(Zorig,'cluster',3,[],method);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, -0.030, -0.369, -0.276, -0.521, 0.243, 0.319, 0.659, -0.517, 1.524, 1.703, 1.286, 1.113; ...
2.900, -0.693, 0.817, -0.209, 0.252, -0.436, 0.675, 0.037, -0.474, 2.714, -1.330, 0.665, -0.767; ...
2.000, -0.157, -0.214, -0.760, -0.153, 0.259, 0.474, -0.637, -0.003, -1.267, -0.478, -1.431, 1.264; ...
3.200, -0.644, 0.260, 0.371, 0.757, 0.589, -0.198, -0.317, 0.652, 0.196, 1.320, -1.912, -1.998];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% Combastel - order 1
Zred = reduce(Zorig,'combastel',1);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, 5.187, 0.000, 0.000, 0.000; ...
2.900, 0.000, 6.362, 0.000, 0.000; ...
2.000, 0.000, 0.000, 4.905, 0.000; ...
3.200, 0.000, 0.000, 0.000, 5.451];
   
% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% Combastel - order 3
Zred = reduce(Zorig,'combastel',3);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, 0.319, -0.276, 0.119, -0.369, -0.517, -0.030, -0.521, 0.659, 2.377, 0.000, 0.000, 0.000; ...
2.900, 0.675, -0.209, -0.646, 0.817, -0.474, -0.693, 0.252, 0.037, 0.000, 2.559, 0.000, 0.000; ...
2.000, 0.474, -0.760, -0.636, -0.214, -0.003, -0.157, -0.153, -0.637, 0.000, 0.000, 1.871, 0.000; ...
3.200, -0.198, 0.371, 0.058, 0.260, 0.652, -0.644, 0.757, -0.317, 0.000, 0.000, 0.000, 2.194];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% constOpt - order 1
w = warning();
warning('off');
Zred = reduce(Zorig_2D,'constOpt',1);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, -2.043, -0.903; ...
2.900, -1.873, 1.611];
   
% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% constOpt - order 3
Zred = reduce(Zorig_2D,'constOpt',3);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, -0.521, 0.319, 0.335, -0.517, -2.043, -0.903; ...
2.900, 0.252, 0.675, 0.340, -0.474, -1.873, 1.611];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);
warning(w);

%% Girard - order 1
Zred = reduce(Zorig,'girard',1);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, 5.187, 0.000, 0.000, 0.000; ...
2.900, 0.000, 6.362, 0.000, 0.000; ...
2.000, 0.000, 0.000, 4.905, 0.000; ...
3.200, 0.000, 0.000, 0.000, 5.451];
   
% check result
res(end+1) = compareMatrices(Zmat,true_mat);


%% Girard - order 3
Zred = reduce(Zorig,'girard',3);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, -0.030, -0.369, -0.276, -0.521, 0.243, 0.319, 0.659, -0.517, 2.253, 0.000, 0.000, 0.000; ...
2.900, -0.693, 0.817, -0.209, 0.252, -0.436, 0.675, 0.037, -0.474, 0.000, 2.769, 0.000, 0.000; ...
2.000, -0.157, -0.214, -0.760, -0.153, 0.259, 0.474, -0.637, -0.003, 0.000, 0.000, 2.248, 0.000; ...
3.200, -0.644, 0.260, 0.371, 0.757, 0.589, -0.198, -0.317, 0.652, 0.000, 0.000, 0.000, 1.663];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% MethA - order 1
Zred = reduce(Zorig,'methA',1);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, 0.912, 1.170, 1.955, 3.604; ...
2.900, -2.489, -2.100, 4.137, 0.202; ...
2.000, 2.779, 1.247, 2.905, -3.483; ...
3.200, -3.608, 2.836, -1.214, -1.733];
   
% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% MethA - order 3
Zred = reduce(Zorig,'methA',3);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, -0.030, -0.369, -0.276, -0.521, 0.243, 0.319, 0.659, -0.517, -1.587, 0.719, -0.200, 0.322; ...
2.900, -0.693, 0.817, -0.209, 0.252, -0.436, 0.675, 0.037, -0.474, -0.057, 2.581, 0.920, -1.750; ...
2.000, -0.157, -0.214, -0.760, -0.153, 0.259, 0.474, -0.637, -0.003, -0.299, -0.667, -0.791, -1.722; ...
3.200, -0.644, 0.260, 0.371, 0.757, 0.589, -0.198, -0.317, 0.652, 1.153, 0.189, -1.448, 0.157];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% MethB - order 1
Zred = reduce(Zorig,'methB',1,6);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, 0.800, -2.483, -3.376, -0.207; ...
2.900, -4.343, 5.497, -3.095, -4.777; ...
2.000, -4.276, -1.440, -0.020, -1.082; ...
3.200, 0.390, 1.749, 4.258, -4.439];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% MethB - order 3
Zred = reduce(Zorig,'methB',3,6);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, -0.030, -0.369, -0.276, -0.521, 0.243, 0.319, 0.659, -0.517, -0.196, 1.687, 0.318, 0.279; ...
2.900, -0.693, 0.817, -0.209, 0.252, -0.436, 0.675, 0.037, -0.474, 0.903, 1.713, -0.738, -1.515; ...
2.000, -0.157, -0.214, -0.760, -0.153, 0.259, 0.474, -0.637, -0.003, -0.777, -1.335, 1.338, -1.492; ...
3.200, -0.644, 0.260, 0.371, 0.757, 0.589, -0.198, -0.317, 0.652, -1.422, 0.121, -1.670, 0.136];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% MethC - order 1
Zred = reduce(Zorig,'methC',1,[12,6]);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, 1.167, 1.972, -3.424, 3.823; ...
2.900, -2.093, 4.173, 1.656, 0.215; ...
2.000, 1.244, 2.930, -1.005, -3.695; ...
3.200, 2.828, -1.224, 4.974, -1.839];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% MethC - order 3
Zred = reduce(Zorig,'methC',3,[12,6]);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, -0.030, -0.369, -0.276, -0.521, 0.243, 0.319, 0.659, -0.517, -0.196, 1.687, 0.318, 0.279; ...
2.900, -0.693, 0.817, -0.209, 0.252, -0.436, 0.675, 0.037, -0.474, 0.903, 1.713, -0.738, -1.515; ...
2.000, -0.157, -0.214, -0.760, -0.153, 0.259, 0.474, -0.637, -0.003, -0.777, -1.335, 1.338, -1.492; ...
3.200, -0.644, 0.260, 0.371, 0.757, 0.589, -0.198, -0.317, 0.652, -1.422, 0.121, -1.670, 0.136];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% PCA - order 1
Zred = reduce(Zorig,'pca',1);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, 0.836, 4.047, -1.271, 2.732; ...
2.900, -6.287, 0.859, -1.147, -0.087; ...
2.000, -1.333, 1.284, 4.633, 0.747; ...
3.200, -0.757, -4.929, -0.037, 2.422];
   
% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% PCA - order 3
Zred = reduce(Zorig,'pca',3);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, -0.030, -0.369, -0.276, -0.521, 0.243, 0.319, 0.659, -0.517, -0.672, -1.903, -1.298, -0.002; ...
2.900, -0.693, 0.817, -0.209, 0.252, -0.436, 0.675, 0.037, -0.474, 2.574, -0.933, 0.021, -0.019; ...
2.000, -0.157, -0.214, -0.760, -0.153, 0.259, 0.474, -0.637, -0.003, 0.616, 1.596, -1.322, 0.606; ...
3.200, -0.644, 0.260, 0.371, 0.757, 0.589, -0.198, -0.317, 0.652, -0.214, -0.650, 0.525, 1.520];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% redistribute - order 1
Zred = reduce(Zorig,'redistribute',1);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);

% true result
true_mat = [ ...
-1.300, -5.551, -0.239, -6.361, 4.192; ...
2.900, -5.090, -5.512, 3.077, 0.235; ...
2.000, -0.032, -1.249, -1.868, -4.052; ...
3.200, 7.001, -5.122, 9.243, -2.016];

% check result
res(end+1) = compareMatrices(Zmat,true_mat);

%% redistribute - order 3
Zred = reduce(Zorig,'redistribute',3);

% obtain zonotope matrix
Zmat = round([Zred.c,Zred.G],3);
% Zmat_sorted = sortMat(Zmat);

% true result
true_mat = [ ...
-1.300, 0.541, 0.400, -0.276, 0.134, -0.517, -0.622, 1.345, -0.369, -0.746, 0.144, 0.313, -0.030; ...
2.900, -0.971, 0.846, -0.209, -0.311, -0.474, 0.301, 0.076, 0.817, 0.256, -0.781, 1.121, -0.693; ...
2.000, 0.577, 0.594, -0.760, 0.564, -0.003, -0.183, -1.300, -0.214, -0.093, -0.769, -0.290, -0.157; ...
3.200, 1.312, -0.248, 0.371, -0.705, 0.652, 0.903, -0.647, 0.260, -0.103, 0.070, 0.082, -0.644];
% true_mat_sorted = sortMat(true_mat);

% check result
% res(end+1) = all(all(Zmat_sorted == true_mat_sorted));
res(end+1) = compareMatrices(Zmat,true_mat);

% plot
% figure; hold on; box on;
% plot(Zorig); plot(Zred,[1,2],'r'); plot(zonotope(true_mat),[1,2],'g');
% close;

res = all(res);

% ------------------------------ END OF CODE ------------------------------
