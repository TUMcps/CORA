function res = test_zonotope_in
% test_zonotope_in - unit test function of in
%
% Syntax:  
%    res = test_zonotope_in
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
% Last update:  14-Sep-2019
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% create parallelotopes
P1 = zonotope([-3.8, -4, 3; 1.2, 3, -4]);
P2 = zonotope([-3.8, -8, 2; 1.2, 10, -10]);

% obtain results
int_1 = in(P1,Z1);
int_2 = in(P2,Z1);

% true results
true_int_1 = 0;
true_int_2 = 1;   

% check result
res = (int_1==true_int_1) & (int_2==true_int_2);

if res
    disp('test_in successful');
else
    disp('test_in failed');
end

%------------- END OF CODE --------------
