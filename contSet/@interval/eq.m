function res = eq(I1,I2,varargin)
% eq - Overloads the '==' operator for exact comparison of two intervals
%
% Syntax:
%    res = I1 == I2
%    res = eq(I1,I2)
%    res = eq(I1,I2,tol)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false whether intervals are equal
%
% Example: 
%    I1 = interval([1; -1], [2; 1]);
%    I2 = interval([1; -2], [2; 2]);
%    I1 == I2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/isequal

% Authors:       Matthias Althoff
% Written:       05-August-2016 
% Last update:   23-December-2022 (MW, call isequal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(I1,I2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
