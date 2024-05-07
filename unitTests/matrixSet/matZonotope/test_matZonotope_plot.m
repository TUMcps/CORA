function resvec = test_matZonotope_plot
% test_matZonotope_plot - unit test function for plot
%
% Syntax:
%    res = test_matZonotope_plot
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
% See also: -

% Authors:       Tobias Ladner
% Written:       26-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

C = [0 2; 1 -1; 1 -2];
G = []; G(:,:,1) = [1 1; -1 0; -2 1]; G(:,:,2) = [-2 0; 0 1; 1 -1];
matZ = matZonotope(C,G);

figure;
plot(matZ)
close;

resvec(end+1) = true;

% combine results
resvec = all(resvec);

% ------------------------------ END OF CODE ------------------------------
