function res = eq(sys1,sys2,varargin)
% eq - overloads '==' operator to check if two linear systems are equal
%
% Syntax:
%    res = sys1 == sys2
%    res = eq(sys1,sys2)
%    res = eq(sys1,sys2,tol)
%
% Inputs:
%    sys1 - linearSys object
%    sys2 - linearSys object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    sys1 = linearSys(1,0);
%    sys2 = linearSys(1,1);
%    res = sys1 == sys2
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

res = isequal(sys1,sys2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
