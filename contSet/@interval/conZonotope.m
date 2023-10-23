function cZ = conZonotope(I)
% conZonotope - Converts an interval object into a conZonotope object
%
% Syntax:
%    cZ = conZonotope(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    I = interval([1;-1], [2;1]);
%    cZ = conZonotope(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Authors:       Niklas Kochdumper
% Written:       21-November-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

cZ = conZonotope(zonotope(I));

% ------------------------------ END OF CODE ------------------------------
