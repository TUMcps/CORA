function pgon = lift_(pgon,N,proj)
% lift_ - lifts a polygon onto a higher-dimensional space
%
% Syntax:
%    pgon = lift_(pgon,N,dims)
%
% Inputs:
%    pgon - polygon object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    pgon - projected polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/lift, contSet/projectHighDim

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if length(proj) == N
    % use project
    pgon = project(pgon,proj);
else
    % projection to higher dimension is not defined as function expects new
    % dimensions to be unbounded
    throw(CORAerror('CORA:notDefined','New dimensions cannot be unbounded as the set representation is always bounded.', 'contSet/projectHighDim'))
end

% ------------------------------ END OF CODE ------------------------------
