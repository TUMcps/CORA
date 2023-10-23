function varargout = size(intMat,varargin)
% size - returns the dimension of the interval matrix
%
% Syntax:
%    n = size(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%    varargin - specified dimension
%
% Outputs:
%    varargout - dimension of the interval matrix
%
% Example: 
%    intMat = intervalMatrix(eye(2),2*eye(2));
%    n = size(intMat)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       17-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
elseif nargin > 1
    [varargout{1:nargout}] = size(intMat.int.inf,varargin{1});
else
    [varargout{1:nargout}] = size(intMat.int.inf);
end

% ------------------------------ END OF CODE ------------------------------
