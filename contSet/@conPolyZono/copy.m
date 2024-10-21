function cPZ_out = copy(cPZ)
% copy - copies the conPolyZono object (used for dynamic dispatch)
%
% Syntax:
%    cPZ_out = copy(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    cPZ_out - copied conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [1 0 1 -1; 0 1 1 1];
%    E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    EC = [0 1 2; 1 0 0; 0 1 0];
%    
%    cPZ = conPolyZono(c,G,E,A,b,EC)
%    cPZ_out = copy(cPZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call copy constructor
cPZ_out = conPolyZono(cPZ);

% ------------------------------ END OF CODE ------------------------------
