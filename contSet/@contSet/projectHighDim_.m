function S = projectHighDim_(S,N,proj,varargin)
% projectHighDim_ - project a set to a higher-dimensional space,
%    having the new dimensions bounded at 0
%
% Syntax:
%    S = projectHighDim_(S,N,proj)
%
% Inputs:
%    S - contSet object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional polytope object
%
% Outputs:
%    S - contSet object in the higher-dimensional space
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/projectHighDim, contSet/project, contSet/lift

% Authors:       Tobias Ladner
% Written:       13-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init identity matrix
I = eye(N);

% project
S = I(:,proj) * S;

% ------------------------------ END OF CODE ------------------------------
