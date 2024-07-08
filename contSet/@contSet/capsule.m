function C = capsule(S,varargin)
% capsule - conversion to capsule objects
%
% Syntax:
%    C = capsule(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    C - capsule object
%
% Example:
%    I = interval([1;2],[4;3]);
%    C = capsule(I);
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
