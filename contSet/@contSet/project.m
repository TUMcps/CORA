function S = project(varargin)
% project - projects a set onto the specified dimensions
%
% Description:
%    computes the set { [s_{(dims_{(1)})}, \dots ,s_{(dims_{(m)})}]^T | s \in \mathcal{S} } \subset \R^m
%
% Syntax:
%    S = project(S, dims)
%
% Inputs:
%    S - contSet object
%    dims - dimensions of projection
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/projectHighDim, contSet/lift

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
