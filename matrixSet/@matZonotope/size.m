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

narginchk(1,2);
if nargin > 1
    [varargout{1:nargout}] = size(matZ.C,varargin{1});
else
    [varargout{1:nargout}] = size(matZ.C);
end

% ------------------------------ END OF CODE ------------------------------
