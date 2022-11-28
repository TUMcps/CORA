function matZ = matZonotope(intMat)
% matZonotope - converts an interval matrix to a matrix zonotope
%
% Syntax:  
%    matZ = matZonotope(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    matZ - zonotope matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      21-June-2010 
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%convert to interval
I=interval(intMat);

%convert to zonotope
Z=zonotope(I);
Z=deleteZeros(Z); %delete zero generators

%convert to matrix zonotope
matZ=matZonotope(Z);

%------------- END OF CODE --------------