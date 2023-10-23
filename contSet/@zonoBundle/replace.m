function zB = replace(zB,index,Z)
% replace - replaces a zonotope at an index position by another zonotope
%
% Syntax:
%    zB = replace(zB,index,Z)
%
% Inputs:
%    zB - zonoBundle object
%    index - index where zonotope is replaced
%    Z - zonotope object
%
% Outputs:
%    zB - zonoBundle object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       01-December-2010
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%replace zonotope
zB.Z{index} = Z;

% ------------------------------ END OF CODE ------------------------------
