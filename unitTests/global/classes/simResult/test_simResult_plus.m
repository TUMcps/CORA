function res = test_simResult_plus
% test_simResult_plus - unit test function for plus
%
% Syntax:
%    res = test_simResult_plus()
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
shift = [0.1;0.2;0.3];
simRes_out = simRes + shift;
resvec(end+1) = isequal(simRes_out.x{1},x{1}+shift') ...
    && isequal(simRes_out.x{2},x{2}+shift');

shift = 2;
simRes_out = simRes + shift;
resvec(end+1) = isequal(simRes_out.x{1},x{1}+shift') ...
    && isequal(simRes_out.x{2},x{2}+shift');

shift = 2;
simRes_out = shift + simRes;
resvec(end+1) = isequal(simRes_out.x{1},x{1}+shift') ...
    && isequal(simRes_out.x{2},x{2}+shift');

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
