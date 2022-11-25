function str = bracketSubs(str)
% bracketSubs - 
%
% Syntax:  
%    str = bracketSubs(str)
%
% Inputs:
%    str - string
%
% Outputs:
%    str - string
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author:       ---
% Written:      ---
% Last update:  01-May-2020 (MW, added header)
% Last revision:---

%------------- BEGIN CODE --------------

%generate left and right brackets
str=strrep(str,'L','(');
str=strrep(str,'R',')');

end

%------------- END OF CODE --------------