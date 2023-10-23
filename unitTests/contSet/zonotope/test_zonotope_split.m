function res = test_zonotope_split
% test_zonotope_split - unit test function of split
%
% Syntax:
%    res = test_zonotope_split
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-14;

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% create halfspace
h = halfspace([1; -1], -2);

% obtain result 1
Zsplit_1 = split(Z1);

% obtain result 2
Zsplit_2 = split(Z1,2);

% obtain result 3
Zsplit_3 = split(Z1,[1; 1]);

% obtain result 4
Zsplit_4 = split(Z1,h);

%SPLIT FUNCTIONS WITH 3 OPERANDS NOT YET TESTED

% obtain zonotope matrix -- split 1
for i=1:length(Zsplit_1)
    for j=1:length(Zsplit_1{1})
        c_1{i}{j} = Zsplit_1{i}{j}.c;
        G_1{i}{j} = Zsplit_1{i}{j}.G;
    end
end

% obtain zonotope matrix -- split 2
for i=1:length(Zsplit_2)
    c_2{i} = Zsplit_2{i}.c;
    G_2{i} = Zsplit_2{i}.G;
end

% obtain zonotope matrix -- split 3
for i=1:length(Zsplit_3)
    c_3{i} = Zsplit_3{i}.c;
    G_3{i} = Zsplit_3{i}.G;
end

% obtain zonotope matrix -- split 4
for i=1:length(Zsplit_4)
    c_4{i} = Zsplit_4{i}.c;
    G_4{i} = Zsplit_4{i}.G;
end

% true result -- split 1
true_c_1{1}{1} = [-7; 1];
true_G_1{1}{1} = [3, 0; 0, 9];
true_c_1{1}{2} = [-1; 1];
true_G_1{1}{2} = [3, 0; 0, 9];
true_c_1{2}{1} = [-4; -3.5];
true_G_1{2}{1} = [6, 0; 0, 4.5];
true_c_1{2}{2} = [-4; 5.5];
true_G_1{2}{2} = [6, 0; 0, 4.5];
        
% true result -- split 2
true_c_2{1} = [-4; -3.5];
true_G_2{1} = [6, 0; 0, 4.5];
true_c_2{2} = [-4; 5.5];
true_G_2{2} = [6, 0; 0, 4.5];
        
% true result -- split 3
true_c_3{1} = [-5.25; -0.25];
true_G_3{1} = [1.25, -2.5, -2.5, -2.5; ...
            1.25, 2.5, 2.5, 2.5];
true_c_3{2} = [-2.75; 2.25];
true_G_3{2} = [-1.25, -2.5, -2.5, -2.5; ...
            -1.25, 2.5, 2.5, 2.5];
        
% true result -- split 4
true_c_4{1} = [-7; 4];
true_G_4{1} = [4.5, -0.5, 0.5, 1.5; ...
            -4.5, -0.5, 0.5, 1.5];
true_c_4{2} = [0.5; -3.5];
true_G_4{2} = [-3, -0.5, 0.5, 1.5; ...
            3, -0.5, 0.5, 1.5];

% check result --split 1
res_1 = true;
for i=1:length(G_1)
    for j=1:length(G_1{1})
        res_1 = res_1 && compareMatrices(true_c_1{i}{j},c_1{i}{j},tol) ...
            && compareMatrices(true_G_1{i}{j},G_1{i}{j},tol);
    end
end

% check result --split 2
res_2 = true;
for i=1:length(G_2)
    res_2 = res_2 && compareMatrices(true_c_2{i},c_2{i},tol) ...
            && compareMatrices(true_G_2{i},G_2{i},tol);
end

% check result --split 3
res_3 = true;
for i=1:length(G_3)
    res_3 = res_3 && compareMatrices(true_c_3{i},c_3{i},tol) ...
            && compareMatrices(true_G_3{i},G_3{i},tol);
end

% check result --split 4
res_4 = true;
for i=1:length(G_4)
    res_4 = res_4 && compareMatrices(true_c_4{i},c_4{i},tol) ...
            && compareMatrices(true_G_4{i},G_4{i},tol);
end

% combined check
res = res_1 & res_2 & res_3 & res_4;

% ------------------------------ END OF CODE ------------------------------
