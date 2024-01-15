function aff = affine(S,varargin)
% affine - conversion to affine objects
%
% Syntax:
%    aff = affine(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    aff - affine object
%
% Example:
%    E = ellipsoid([1 0; 0 2];[1;-1]);
%    aff = affine(E);
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
