function [c] = center(qZ)
% center - Returns the center of a qudZonotope
%
% Syntax:  
%    [c] = center(qZ)
%
% Inputs:
%    qZ - quadZonotope object
%
% Outputs:
%    c - center of the quadZonotope qZ
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
% Written:      05-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

c=qZ.c;

%------------- END OF CODE --------------