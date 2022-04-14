function [Zbundle] = replace(Zbundle,index,Z)
% replace - replaces a zonotope at an index position by another zonotope
%
% Syntax:  
%    [Zbundle] = replace(Zbundle,index,Z)
%
% Inputs:
%    Zbundle - zonotope bundle
%    index - index where zonotope is replaced
%    Z - zonotope
%
% Outputs:
%    Zbundle - zonotope bundle
%
% Example: 
%    ---
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      01-December-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%replace zonotope
Zbundle.Z{index}=Z;


%------------- END OF CODE --------------