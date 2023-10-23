function vol = volume_(varargin)
% volume_ - computes the volume of a set
%    (internal use, see also contSet/volume)
%
% Syntax:
%    V = volume_(S,method)
%
% Inputs:
%    S - contSet object
%    method - (optional) method
%
% Outputs:
%    V - numeric matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/volume

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
