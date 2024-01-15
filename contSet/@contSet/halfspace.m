function hs = halfspace(S,varargin)
% halfspace - conversion to halfspace objects
%
% Syntax:
%    hs = halfspace(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    hs - halfspace object
%
% Example:
%    E = ellipsoid([1 0; 0 2];[1;-1]);
%    hs = halfspace(E);
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
