function zB = zonoBundle(S,varargin)
% zonoBundle - conversion to zonoBundle objects
%
% Syntax:
%    zB = zonoBundle(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    zB - zonoBundle object
%
% Example:
%    Z = zonotope([1;1],[1 1 1; 1 -1 0]);
%    zB = zonoBundle(Z);
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
