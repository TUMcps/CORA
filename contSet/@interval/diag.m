function res = diag(I, varargin)
% diag - Create diagonal matrix or get diagonal elements of matrix
%
% Syntax:
%    res = diag(I)
%    res = diag(I,k)
%
% Inputs:
%    I - interval object
%    k - (optional) diagonal number
%
% Outputs:
%    res - diagonal matrix or diagonal elements of matrix
%
% Example: 
%    I = interval([1; -1], [2; 0]);
%    d = diag(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: diag

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       02-November-2017 
% Last update:   18-July-2023 (TL, getter for k-th diagonal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = dim(I);
if numel(n) > 2
    throw(CORAerror('CORA:wrongValue','first','Interval must not be an n-d array with n > 2.'))
end

% obtain result
res = interval(diag(I.inf, varargin{:}), diag(I.sup, varargin{:}));

% ------------------------------ END OF CODE ------------------------------
