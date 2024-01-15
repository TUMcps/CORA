function ls = levelSet(S,varargin)
% levelSet - conversion to levelSet objects
%
% Syntax:
%    ls = levelSet(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    ls - levelSet object
%
% Example:
%    E = ellipsoid([1 0; 0 2];[1;-1]);
%    ls = levelSet(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       23-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror('CORA:noops',S));

% ------------------------------ END OF CODE ------------------------------
