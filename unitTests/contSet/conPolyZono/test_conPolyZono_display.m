function res = test_conPolyZono_display
% test_conPolyZono_display - unit test function of display
%
% Syntax:
%    res = test_conPolyZono_display
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
% Written:       04-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% test empty set
cPZ = conPolyZono.empty(3);
display(cPZ);
res(end+1) = true;

% test example in doctring
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];

cPZ = conPolyZono(c,G,E,A,b,EC);

% check display function
display(cPZ);
res(end+1) = true;

% check direct output
cPZ
res(end+1) = true;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
