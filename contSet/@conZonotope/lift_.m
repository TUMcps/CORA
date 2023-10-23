function cZ = lift_(cZ,N,proj)
% lift_ - lifts a conZonotope onto a higher-dimensional space
%
% Syntax:
%    cZ = lift_(cZ,N,dims)
%
% Inputs:
%    cZ - conZonotope object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    cZ - projected conZonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/lift, contSet/projectHighDim

% Authors:       Tobias Ladner
% Written:       19-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if length(proj) == N
    % use project
    cZ = project(cZ,proj);
else
    % projection to higher dimension is not defined as function expects new
    % dimensions to be unbounded
    throw(CORAerror('CORA:notDefined','New dimensions cannot be unbounded as the set representation is always bounded.', 'contSet/projectHighDim'))
end

% ------------------------------ END OF CODE ------------------------------
