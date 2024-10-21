function res = eq(nlnsys1,nlnsys2,varargin)
% eq - overloads '==' operator to check if two nonlinear systems are equal
%
% Syntax:
%    res = nlnsys1 == nlnsys2
%    res = eq(nlnsys1,nlnsys2)
%    res = eq(nlnsys1,nlnsys2,tol)
%
% Inputs:
%    nlnsys1 - nonlinearSys object
%    nlnsys2 - nonlinearSys object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example:
%    f = @(x,u) [x(1)^2 - u(1); x(2)];
%    g = @(x,u) [x(2)^2 - u(1); x(1)]; 
%    nlnsys1 = nonlinearSys(f);
%    nlnsys2 = nonlinearSys(g);
%    res = nlnsys1 == nlnsys2
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

res = isequal(nlnsys1,nlnsys2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
