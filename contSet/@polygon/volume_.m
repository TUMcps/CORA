function vol = volume_(pgon,varargin)
% volume_ - computes the volume of a polygon
%
% Syntax:
%    vol = volume_(pgon,method,order)
%
% Inputs:
%    pgon - polygon object
%    method - (optional) method for approximation 
%    order - (optional) zonotope order for reduction before computation
%
% Outputs:
%    vol - volume
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/volume

% Authors:       Tobias Ladner
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

vol = pgon.set.area;

end

% ------------------------------ END OF CODE ------------------------------
