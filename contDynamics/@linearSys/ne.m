function res = ne(linsys1,linsys2,varargin)
% ne - overloads '~=' operator to check if two linear systems are not equal
%
% Syntax:
%    res = linsys1 ~= linsys2
%    res = ne(linsys1,linsys2)
%    res = ne(linsys1,linsys2,tol)
%
% Inputs:
%    linsys1 - linearSys object
%    linsys2 - linearSys object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    linsys1 = linearSys(1,0);
%    linsys2 = linearSys(1,1);
%    res = linsys1 ~= linsys2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isequal(linsys1,linsys2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
