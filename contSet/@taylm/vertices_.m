function V = vertices_(tay,varargin)
% vertices_ - returns the vertices of a Taylor model
%
% Syntax:
%    V = vertices_(tay)
%
% Inputs:
%    obj - a taylm object
%
% Outputs:
%    res - a taylm object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices, polyZonotope/vertices_

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to polyZonotope
pZ = polyZonotope(tay);

% compute vertices
V = vertices(pZ,varargin{:});

end

% ------------------------------ END OF CODE ------------------------------
