function E = ellipsoid(S,varargin)
% ellipsoid - conversion to ellipsoid objects
%
% Syntax:
%    E = ellipsoid(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    E - ellipsoid object
%
% Example:
%    pZ = polyZonotope([0;1],[2 0 1;0 2 1],[0;0.5],[1 0 3;0 1 1]);
%    E = ellipsoid(pZ);
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
