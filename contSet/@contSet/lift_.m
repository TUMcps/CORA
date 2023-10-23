function S = lift_(varargin)
% lift_ - lifts a set to a higher-dimensional space,
%    having the new dimensions unbounded
%
% Syntax:
%    S = lift_(S,N,proj)
%
% Inputs:
%    S - contSet object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional contSet object
%
% Outputs:
%    S - contSet object in the higher-dimensional space
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/lift, contSet/project, contSet/projectHighDim

% Authors:       Tobias Ladner
% Written:       19-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
