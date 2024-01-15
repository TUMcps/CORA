function I = interval(S,varargin)
% interval - conversion to interval objects
%
% Syntax:
%    I = interval(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    I - interval object
%
% Example:
%    E = ellipsoid([1 0; 0 2];[1;-1]);
%    I = interval(E);
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
