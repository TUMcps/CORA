function res = test_zonotope_enclose
% test_zonotope_enclose - unit test function of enclose
%
% Syntax:  
%    res = test_zonotope_enclose
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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020 (MW, enhance randomness)
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Analytical Test ------------------------------------------------------

% create zonotopes
Z1 = zonotope([1,2,3,4; 5 6 7 8]);
Z2 = zonotope([9, 10, 11; 12, 13, 14]);

% obtain enclosing zonotope
Z3 = enclose(Z1,Z2);

% obtain zonotope matrix
Zmat = Z3.Z;

% true result
true_mat = [5, 6, 7, -4, -4, -4, 4; ...
            8.5, 9.5, 10.5, -3.5, -3.5, -3.5, 8];

% check result
res_val = all(all(Zmat == true_mat));


% add results
res = res_val;

if res
    disp('test_enclose successful');
else
    disp('test_enclose failed');
end

%------------- END OF CODE --------------
