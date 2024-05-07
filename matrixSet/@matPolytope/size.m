function varargout = size(matP,varargin)
% size - returns the dimension of the matrix polytope
%
% Syntax:
%    n = size(matP)
%
% Inputs:
%    matP - matPolytope object
%    varargin - specified dimension
%
% Outputs:
%    varargout - dimension of the matrix polytope
%
% Example:
%    matP = matPolytope({[1 2; 0 1],[1 3; -1 2]});
%    n = size(matP)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       17-January-2023
% Last update:   02-May-2024 (TL, simplified)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin == 1
    dim = 1:2;
elseif nargin == 2
    dim = varargin{1};
else
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% read size of vertices
[varargout{1:nargout}] = size(matP.V,dim);

end

% ------------------------------ END OF CODE ------------------------------
