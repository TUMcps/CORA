function hyp = conHyperplane(S,varargin)
% conHyperplane - conversion to conHyperplane objects
%
% Syntax:
%    hyp = conHyperplane(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    hyp - conHyperplane object
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
