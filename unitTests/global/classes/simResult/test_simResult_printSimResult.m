function res = test_simResult_printSimResult
% test_simResult_printSimResult - unit test function of printSpec
%
% Syntax:
%    res = test_simResult_printSimResult
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
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test print of a simple system
N = 10; n = 3;
x = {rand(N,n),rand(N,n)};
t = {(1:N)',(1:N)'};
simRes = simResult(x,t);

printSimResult(simRes);
printSimResult(simRes,'high');
printSimResult(simRes,'high',true);
printSimResult(simRes,'high',false);

% test fid
filename = 'test.txt';
printSimResult(filename,simRes,'high',true);
simRes_copy = eval(fileread(filename));
% assert(isequal(simRes,simRes_copy)); % not implemented
delete(filename)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
