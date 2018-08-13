function [qZ] = project(qZ,dim)
% project - Returns a quadZonotope which is projected onto the specified
% dimensions
%
% Syntax:  
%    [qZ] = project(qZ,dim)
%
% Inputs:
%    qZ - quadZonotope object
%    dim - projected dimensions
%
% Outputs:
%    qZ - quadZonotope
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
% Written:      04-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

qZ.c=qZ.c(dim,:);
qZ.G=qZ.G(dim,:);
qZ.Gsquare=qZ.Gsquare(dim,:);
qZ.Gquad=qZ.Gquad(dim,:);
qZ.Grest=qZ.Grest(dim,:);

%------------- END OF CODE --------------