function res = test_zonotope_interval
% test_zonotope_interval - unit test function of interval
%
% Syntax:  
%    res = test_zonotope_interval
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
I1 = interval(Z1);

% obtain results
inf = infimum(I1);
sup = supremum(I1);

% true results
true_inf = [-10; -8];
true_sup = [2; 10];   

% check result
res = all(inf==true_inf) & all(sup==true_sup);

if res
    disp('test_zonotope_interval successful');
else
    disp('test_zonotope_interval failed');
end

%------------- END OF CODE --------------
