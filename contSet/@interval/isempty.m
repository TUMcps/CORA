function res = isempty(I)
% isempty - checks if an interval is the empty set
%
% Syntax:
%    res = isempty(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - true/false
%
% Example: 
%    I = interval([-1;-2],[3;4]);
%    res = isempty(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      12-December-2010
% Last update:  19-November-2015 (Daniel Althoff)
% Last revision:---

%------------- BEGIN CODE --------------

%return result
res = isempty(I.inf) || isempty(I.sup);

%------------- END OF CODE --------------