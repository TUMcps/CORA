function res = test_zonotope_quadMap
% test_zonotope_quadMap - unit test function of quadMap
%
% Syntax:  
%    res = test_zonotope_quadMap
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

% create zonotopes
Z1 = zonotope([-4, -3, -2; 1, 2, 3]);
Z2 = zonotope([1, 4, 2; -3, 2, -1]);

% create matrices
Q{1} = [1 -2; -3 4];
Q{2} = [0.5 0; 2 -1];

% 1. quadMapSingle: Z1*Q*Z1 -----------------------------------------------

% obtain result
Zres = quadMap(Z1,Q);

% obtain zonotope matrix
Zmat = Zres.Z;

% true result
true_mat = [102.5, 27.5, 35, 95, 110, 125; ...
            -16.25, -5.75, -9.5, -14, -26, -32];

% compare solutions
res(1) = all(all(Zmat == true_mat));

% 2. quadMapMixed: Z1*Q*Z2 ------------------------------------------------

% obtain result
Zres = quadMap(Z1,Z2,Q);

% obtain zonotope matrix
Zmat = Zres.Z;

% true result
true_mat = [-43, -51, -59, -4, -8, -12, -26, -32, -38;
            3, 8.5, 14, -2, 6, 14, 1, 7, 13];

% compare solutions
res(2) = all(all(Zmat == true_mat));


% gather results
res = all(res);

if res
    disp('test_zonotope_quadMap successful');
else
    disp('test_zonotope_quadMap failed');
end

%------------- END OF CODE --------------
