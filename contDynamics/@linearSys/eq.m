function res = eq(linsys1,linsys2,varargin)
% eq - overloads '==' operator to check if two linear systems are equal
%
% Syntax:
%    res = linsys1 == linsys2
%    res = eq(linsys1,linsys2)
%    res = eq(linsys1,linsys2,tol)
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
%    res = linsys1 == linsys2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-January-2023
% Last update:   10-January-2023 (MW, move code to isequal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(linsys1,linsys2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
