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

% Authors:       Mark Wetzlinger
% Written:       17-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% take first vertex (all vertices have to have same dimension)
if nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
elseif nargin > 1
    if isempty(matP)
        [varargout{1:nargout}] = 0;
    else
        [varargout{1:nargout}] = size(matP.vertex{1},varargin{1});
    end
else
    if isempty(matP)
        [varargout{1:nargout}] = [0,0];
    else
        [varargout{1:nargout}] = size(matP.vertex{1});
    end
end

% ------------------------------ END OF CODE ------------------------------
