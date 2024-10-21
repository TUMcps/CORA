function V = vertices_(pgon,varargin)
% vertices_ - returns the vertices of the polygon
%
% Syntax:
%    V = vertices_(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    V - numeric
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read vertices
V = pgon.set.Vertices';

% sort vertices of multiple regions or holes (separated by nan values)
% V = aux_sortMultipleRegionsAndHoles(V);

% polygons are stored with buffers as the underlying polyshape 
% cannot deal with boundary-only sets, i.e., a point or a line.
% Merging buffered points here again..

if~isempty(V) && ~any(isnan(V),"all")
    % gather close subsequent points 
    V_rounded = unique(round(V',6),'rows')';

    % 
    if size(V_rounded,2) <= 2 % max two points left (thus, point or line)
        % assume that points were buffered in constructor
        V = V_rounded;
    end
end

end


% Auxiliary functions -----------------------------------------------------


% ------------------------------ END OF CODE ------------------------------
