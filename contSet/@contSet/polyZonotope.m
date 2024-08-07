function pZ = polyZonotope(S,varargin)
% polyZonotope - conversion to polyZonotope objects
%
% Syntax:
%    pZ = polyZonotope(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example:
%    E = ellipsoid([1 0; 0 2],[1;-1]);
%    pZ = polyZonotope(E);
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
