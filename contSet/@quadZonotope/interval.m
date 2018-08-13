function [IH] = interval(qZ)
% interval - Overapproximates a quadZonotope by an interval
%
% Syntax:  
%    [IH] = interval(qZ)
%
% Inputs:
%    qZ - quadZonotope object
%
% Outputs:
%    IH - interval object
%
% Example: 
%    ---
%
% Other m-files required: intervalhull(constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      05-September-2012
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%instantiate interval 
IH=interval(zonotope(qZ));

%------------- END OF CODE --------------