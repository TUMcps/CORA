function res = ne(I1,I2,varargin)
% ne - overloads '~='-operator for intervals
%
% Syntax:
%    res = I1 ~= I2
%    res = ne(I1,I2)
%    res = ne(I1,I2,tol)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    I1 = interval(-1,0);
%    I2 = interval(-2,-1);
%    I1 ~= I2
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/isequal

% Authors:       Dmitry Grebenyuk
% Written:       06-August-2017
% Last update:   23-December-2022 (MW, call isequal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isequal(I1,I2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
