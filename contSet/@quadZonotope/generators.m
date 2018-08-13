function [G,Gsquare,Gquad,Grest] = generators(qZ)
% generators - Returns the generators of a quadZonotope
%
% Syntax:  
%    [G,Gsquare,Gquad,Grest] = generators(qZ)
%
% Inputs:
%    qZ - quadZonotope object
%
% Outputs:
%    G,Gsquare,Gquad,Grest - matrices of generators
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      06-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

G=qZ.G;
Gsquare=qZ.Gsquare;
Gquad=qZ.Gquad;
Grest=qZ.Grest;

%------------- END OF CODE --------------