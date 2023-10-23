function res = test_polyZonotope_restructure
% test_polyZonotope_restructure - unit test function for over-approximative
%    polynomial zonotope restructuring
%
% Syntax:
%    res = test_polyZonotope_restructure
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

% Authors:       Niklas Kochdumper
% Written:       29-March-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% create polynomial zonotope
% pZ = polyZonotope([0;0],[0 4 1 -1; 1 2 -1 -1],[-7 1 1;15 1 -1],[1 0 0 1;0 1 3 2]);
pZ = polyZonotope([0;0],[0 4 1 -1 2; 1 2 -1 -1 1],[-7 1 1;15 1 -1],[1 0 0 0 1;0 1 0 3 2; 0 0 1 1 0]);

% restructure the polynomial zonotope
pZres = restructure(pZ,'reduceGirard',2);

% define ground truth
c = [0;0];
G = [4 1 11 0 -1; 2 -1 0 19 -1];
E = [1 0 0 0 3; 0 1 0 0 1; 0 0 1 0 0; 0 0 0 1 0];

% check for correctness
if ~all(withinTol(c,pZres.c)) || ~compareMatrices([G;E],[pZres.G;pZres.E]) ...
        || ~isempty(pZres.GI)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
