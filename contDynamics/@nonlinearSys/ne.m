function res = ne(sys1,sys2,varargin)
% ne - overloads '~=' operator to check if two nonlinear systems are not
%   equal
%
% Syntax:
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
%    f = @(x,u) [x(1)^2 - u(1); x(2)];
%    g = @(x,u) [x(2)^2 - u(1); x(1)]; 
%    sys1 = nonlinearSys(f);
%    sys2 = nonlinearSys(g);
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
