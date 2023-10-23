function cPZ = lift_(cPZ,N,proj)
% lift_ - lifts a conPolyZono onto a higher-dimensional space
%
% Syntax:
%    cPZ = lift_(cPZ,N,dims)
%
% Inputs:
%    cPZ - conPolyZono object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    C - projected conPolyZono
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
    cPZ = project(cPZ,proj);
else
    % projection to higher dimension is not defined as function expects new
    % dimensions to be unbounded
    throw(CORAerror('CORA:notDefined','New dimensions cannot be unbounded as the set representation is always bounded.', 'contSet/projectHighDim'))
end

% ------------------------------ END OF CODE ------------------------------
