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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Niklas Kochdumper
% Written:      29-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% ANALYTICAL TESTS

% TEST 1

% create polynomial zonotope
% pZ = polyZonotope([0;0],[0 4 1 -1; 1 2 -1 -1],[-7 1 1;15 1 -1],[1 0 0 1;0 1 3 2]);
pZ = polyZonotope([0;0],[0 4 1 -1 2; 1 2 -1 -1 1],[-7 1 1;15 1 -1],[1 0 0 0 1;0 1 0 3 2; 0 0 1 1 0]);

% restructure the polynomial zonotope
pZres = restructure(pZ,'reduceGirard',2);

% define ground truth
c = [0;0];
G = [4 1 11 0 -1; 2 -1 0 19 -1];
expMat = [1 0 0 0 3; 0 1 0 0 1; 0 0 1 0 0; 0 0 0 1 0];

% check for correctness
if any(c-pZres.c) || any(any(G-pZres.G)) || ...
   any(any(expMat-pZres.expMat)) || ~isempty(pZres.Grest)
    error('test_polyZonotope_restructure: analytical test 1 failed!');
end





res = true;

%------------- END OF CODE --------------