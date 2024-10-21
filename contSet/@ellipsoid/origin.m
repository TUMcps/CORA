function E = origin(n)
% origin - instantiates an ellipsoid that contains only the origin
%
% Syntax:
%    E = ellipsoid.origin(n)
%
% Inputs:
%    n - dimension (integer, >= 1)
%
% Outputs:
%    E - ellipsoid representing the origin
%
% Example: 
%    E = ellipsoid.origin(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       21-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({{n,'att','numeric',{'scalar','positive','integer'}}});
E = ellipsoid(zeros(n,n),zeros(n,1));

% ------------------------------ END OF CODE ------------------------------
