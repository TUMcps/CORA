function res = test_zonotope_vertices
% test_zonotope_vertices - unit test function of vertices
%
% Syntax:  
%    res = test_zonotope_vertices
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
% Last update:  09-February-2017
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% obtain result
V1 = vertices(Z1);

% obtain matzrix of vertices
Vmat = get(V1,'V');

% true result
true_Vmat = [-8, 0, 2, -4, -4, -10; ...
   2, 0, -8, -4, 6, 10];   

% check results; order of vectors does not matter
pointExists = 1;
i = 1;
while pointExists && i<=length(Vmat(1,:))
    pointExists = 0;
    for j = 1:length(true_Vmat(1,:))
        if all(abs(Vmat(:,i)-true_Vmat(:,j)) < 1e-13)
            pointExists = 1;
        end
    end
    % increment counter i
    i = i+1;
end

% res
res = pointExists;


if res
    disp('test_zonotope_vertices successful');
else
    disp('test_zonotope_vertices failed');
end

%------------- END OF CODE --------------
