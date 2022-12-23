function res = test_zonotope_deleteZeros
% test_zonotope_deleteZeros - unit test function of deleteZeros
%    this encompasses checking the function nonzeroFilter
%
% Syntax:  
%    res = test_zonotope_deleteZeros
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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020 (MW, enhance randomness)
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
c = [1;5];
G = [2,0,4; 6 0 0];
Z = zonotope(c,G);

% obtain zonotope without zeros
Z_ = deleteZeros(Z);
G_ = generators(Z_);

% true result
true_mat = [2, 4; 6, 0];

% check result
res = compareMatrices(G_,true_mat);

%------------- END OF CODE --------------
