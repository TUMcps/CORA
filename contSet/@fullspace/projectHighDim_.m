function fs = projectHighDim_(fs,N,proj)
% projectHighDim_ - projects a full-dimensional space onto a higher-dimensional space
%    case R^0: undefined
%
% Syntax:
%    fs = projectHighDim_(fs,N,dims)
%
% Inputs:
%    fs - fullspace object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    fs - projected fullspace
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/projectHighDim, contSet/lift

% Authors:       Tobias Ladner
% Written:       19-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if length(proj) == N
    % keep as is
else
    % projection to higher dimension is not defined as function expects new
    % dimensions to be bounded at 0
    throw(CORAerror('CORA:notDefined','Cannot bound new dimensions at 0', 'contSet/lift'))
end

% ------------------------------ END OF CODE ------------------------------
