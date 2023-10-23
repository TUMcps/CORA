function c = center(varargin)
% center - returns the center of a set
%
% Syntax:
%    S = center(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    c - numeric
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
