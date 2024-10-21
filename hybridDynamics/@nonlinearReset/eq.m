function res = eq(nonlinReset1,nonlinReset2,varargin)
% eq - checks if two nonlinear reset functions, or a nonlinear reset
%    function and a linear reset function, are equal up to some tolerance
%
% Syntax:
%    res = nonlinReset1 == nonlinReset2
%    res = eq(nonlinReset1,nonlinReset2)
%    res = eq(nonlinReset1,nonlinReset2,tol)
%
% Inputs:
%    nonlinReset1 - nonlinearReset object, linearReset object
%    nonlinReset2 - nonlinearReset object, linearReset object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    f = @(x,u) [-x(1)*x(2); x(2) - u(1)];
%    g = @(x,u) [-x(1)*x(2); x(1) - u(1)];
%    nonlinReset1 = nonlinearReset(f);
%    nonlinReset2 = nonlinearReset(g);
%    nonlinReset1 == nonlinReset1;
%    nonlinReset1 == nonlinReset2;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearReset/isequal

% Authors:       Mark Wetzlinger
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% redirect to isequal
res = isequal(nonlinReset1,nonlinReset2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
