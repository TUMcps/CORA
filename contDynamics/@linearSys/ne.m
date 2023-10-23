function res = ne(sys1,sys2,varargin)
% ne - overloads '~=' operator to check if two linear systems are not equal
%
% Syntax:
%    res = sys1 ~= sys2
%    res = ne(sys1,sys2)
%    res = ne(sys1,sys2,tol)
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
%    res = sys1 ~= sys2
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

res = ~isequal(sys1,sys2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
