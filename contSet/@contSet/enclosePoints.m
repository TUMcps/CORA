function S = enclosePoints(varargin)
% enclosePoints - enclose a point cloud with a set
%
% Syntax:
%    S = contSet.enclosePoints(points)
%    S = contSet.enclosePoints(points, method)
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%    method - (optional) method
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
