function varargout = size(matZ,varargin)
% size - returns the dimension of the matrix zonotope
%
% Syntax:
%    n = size(matZ)
%
% Inputs:
%    matZ - matZonotope object
%    varargin - specified dimension
%
% Outputs:
%    varargout - dimension of the matrix zonotope
%
% Example:
%    C = [0 0; 0 0];
%    G{1} = [1 3; -1 2];
%    G{2} = [2 0; 1 -1];
%
%    matZ = matZonotope(C,G);
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
    [varargout{1:nargout}] = size(matZ.center,varargin{1});
else
    [varargout{1:nargout}] = size(matZ.center);
end

% ------------------------------ END OF CODE ------------------------------
