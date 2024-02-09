function res = test_simResult_mtimes
% test_simResult_mtimes - unit test function for mtimes
%
% Syntax:
%    res = test_simResult_mtimes()
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

% Authors:       Tobias Ladner
% Written:       02-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% init simResult object
n = 3;
N = 100;

x = {rand(N,n),rand(N,n)};
t = {(1:N)',(1:N)'};

simRes = simResult(x,t);

% simple cases
A = [2 3 -1; 1 2 4];
simRes_out = A * simRes;
resvec(end+1) = isequal(simRes_out.x{1},x{1}*A') ...
    && isequal(simRes_out.x{2},x{2}*A');

A = 3;
simRes_out = A * simRes;
resvec(end+1) = isequal(simRes_out.x{1},x{1}*A') ...
    && isequal(simRes_out.x{2},x{2}*A');

A = 3;
simRes_out = simRes * A;
resvec(end+1) = isequal(simRes_out.x{1},x{1}*A') ...
    && isequal(simRes_out.x{2},x{2}*A');

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
