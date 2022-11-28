function c = center(I)
% center - returns the center of an interval
%
% Syntax:  
%    c = center(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    c - center of interval (vector)
%
% Example: 
%    I = interval([-1;1],[1;2]);
%    c = center(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      26-June-2015
% Last update:  02-Sep-2019 (rename mid -> center)
% Last revision:---

%------------- BEGIN CODE --------------
if isempty(I)
    c = []; return ;
end

c = 0.5*(I.inf + I.sup);

%------------- END OF CODE --------------