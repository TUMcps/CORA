function res = test_simResult_monitorSTL
% test_simResult_monitorSTL - unit test function for verifySTL
%
% Syntax:
%    res = test_simResult_monitorSTL()
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
% See also: none

% Authors:       Benedikt Seidl
% Written:       14-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% simulation result
x = [1.1 1; 2.2 2; 3.3 3; 4.4 4; 5.5 5; 6.6 6; 7.7 7; 8.8 8; 9.9 9];
t = [0; 1; 2; 3; 4; 5; 6; 7; 8];

simRes = simResult({x}, {t});

% STL formulas
x = stl('x', 2);

phi{1} = globally(x(1) > x(2), interval(0,8));
phi{2} = finally(x(1) > 5, interval(4,7));

res = true;

for i=1:length(phi)
    res = res && monitorSTL(simRes,phi{i});
end

% ------------------------------ END OF CODE ------------------------------
