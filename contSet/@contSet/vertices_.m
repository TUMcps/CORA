function V = vertices_(varargin)
% vertices_ - computes the vertices of a set
%    (internal use, see also contSet/vertices)
%
% Syntax:
%    V = vertices_(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    V - numeric matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
