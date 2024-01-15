function probZ = probZonotope(S,varargin)
% probZonotope - conversion to probZonotope objects
%
% Syntax:
%    probZ = probZonotope(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    probZ - probZonotope object
%
% Example:
%    E = ellipsoid([1 0; 0 2];[1;-1]);
%    probZ = probZonotope(E);
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
