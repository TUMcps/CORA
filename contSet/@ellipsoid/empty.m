function E = empty(n)
% empty - instantiates an empty ellipsoid
%
% Syntax:
%    E = ellipsoid.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    E - empty ellipsoid
%
% Example: 
%    E = ellipsoid.empty(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if n <= 0
    throw(CORAerror('CORA:wrongValue','first','positive'));
end

% note: we cannot instantiate an n-dimensional empty matrix, thus zeros(0,0)
E = ellipsoid(zeros(0,0),zeros(n,0));

% ------------------------------ END OF CODE ------------------------------
