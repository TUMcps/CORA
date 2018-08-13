function res = test_zonotope_inParallelotope
% test_inParallelotope - unit test function of inParallelotope
%
% Syntax:  
%    res = test_zonotope_inParallelotope
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

% create parallelotopes
P1 = zonotope([-3.8, -4, 3; 1.2, 3, -4]);
P2 = zonotope([-3.8, -8, 2; 1.2, 10, -10]);

% obtain results
int_1 = inParallelotope(Z1,P1);
int_2 = inParallelotope(Z1,P2);

% true results
true_int_1 = 0;
true_int_2 = 1;   

% check result
res = (int_1==true_int_1) & (int_2==true_int_2);

if res
    disp('test_inParallelotope successful');
else
    disp('test_inParallelotope failed');
end

%------------- END OF CODE --------------
