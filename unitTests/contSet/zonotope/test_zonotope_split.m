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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff
% Written:      26-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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
        Zmat_1{i}{j} = get(Zsplit_1{i}{j},'Z');
    end
end

% obtain zonotope matrix -- split 2
for i=1:length(Zsplit_2)
    Zmat_2{i} = get(Zsplit_2{i},'Z');
end

% obtain zonotope matrix -- split 3
for i=1:length(Zsplit_3)
    Zmat_3{i} = get(Zsplit_3{i},'Z');
end

% obtain zonotope matrix -- split 4
for i=1:length(Zsplit_4)
    Zmat_4{i} = get(Zsplit_4{i},'Z');
end

% true result -- split 1
true_mat_1{1}{1} = [-7, 3, 0; ...
            1, 0, 9];
true_mat_1{1}{2} = [-1, 3, 0; ...
            1, 0, 9];
true_mat_1{2}{1} = [-4, 6, 0; ...
            -3.5, 0, 4.5];
true_mat_1{2}{2} = [-4, 6, 0; ...
            5.5, 0, 4.5];
        
% true result -- split 2
true_mat_2{1} = [-4, 6, 0; ...
            -3.5, 0, 4.5];
true_mat_2{2} = [-4, 6, 0; ...
            5.5, 0, 4.5];
        
% true result -- split 3
true_mat_3{1} = [-5.25, 1.25, -2.5, -2.5, -2.5; ...
            -0.25, 1.25, 2.5, 2.5, 2.5];
true_mat_3{2} = [-2.75, -1.25, -2.5, -2.5, -2.5; ...
            2.25, -1.25, 2.5, 2.5, 2.5];
        
% true result -- split 4
true_mat_4{1} = [-7, 4.5, -0.5, 0.5, 1.5; ...
            4, -4.5, -0.5, 0.5, 1.5];
true_mat_4{2} = [0.5, -3, -0.5, 0.5, 1.5; ...
            -3.5, 3, -0.5, 0.5, 1.5];

% check result --split 1
res_1 = 1;
for i=1:length(Zmat_1)
    for j=1:length(Zmat_1{1})
        res_1 = res_1 & all(all(Zmat_1{i}{j} == true_mat_1{i}{j}));
    end
end

% check result --split 2
res_2 = 1;
for i=1:length(Zmat_2)
    res_2 = res_2 & all(all(Zmat_2{i} == true_mat_2{i}));
end

% check result --split 3
res_3 = 1;
for i=1:length(Zmat_3)
    res_3 = res_3 & all(all(abs(Zmat_3{i}-true_mat_3{i}) < 1e-15));
end

% check result --split 4
res_4 = 1;
for i=1:length(Zmat_4)
    res_4 = res_4 & all(all(abs(Zmat_4{i}-true_mat_4{i}) < 1e-15));
end

% combined check
res = res_1 & res_2 & res_3 & res_4;


if res
    disp('test_zonotope_reduce successful');
else
    disp('test_zonotope_reduce failed');
end

%------------- END OF CODE --------------
