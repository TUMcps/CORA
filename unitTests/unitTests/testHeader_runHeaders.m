function res = testHeader_runHeaders
% testHeader_runHeaders - checks if all headers run successfully
%
% Syntax:
%    res = testHeader_runHeaders
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false 
%

% Authors:       Tobias Ladner
% Written:       04-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp(' ')
assert(runHeaders(true));
disp(' ')

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
