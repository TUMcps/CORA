function res = test_simResult_uminus
% test_simResult_uminus - unit test function for uminus
%
% Syntax:
%    res = test_simResult_uminus()
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
simRes_out = -simRes;
resvec(end+1) = isequal(simRes_out.x{1},-x{1}) ...
    && isequal(simRes_out.x{2},-x{2});

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
