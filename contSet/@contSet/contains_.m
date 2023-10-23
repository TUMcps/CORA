function res = contains_(varargin)
% contains_ - determines if a set contains another set or a point
%    (internal use, see also contSet/contains)
%
% Syntax:
%    res = contains_(S1, S2)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object or numeric
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
