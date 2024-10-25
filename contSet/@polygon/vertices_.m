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

% check if vertices are stored directly
if ~isempty(pgon.V)
    V = pgon.V;
    return
end

% read vertices
V = pgon.set.Vertices';

% sort vertices of multiple regions or holes (separated by nan values)
% V = aux_sortMultipleRegionsAndHoles(V);

% polygons are stored with buffers as the underlying polyshape 
% cannot deal with boundary-only sets, i.e., a point or a line.
% Merging buffered points here again..

if~isempty(V) && ~any(isnan(V),"all") && size(V,2) < 10
    % assume that points were buffered in constructor (points/line)
    V_buffered = V;
    V_true = [];
    while ~isempty(V_buffered) && size(V_buffered,2) < 10
        % find vertices close to the (current) first vertex (+/- TOL + eps)
        idx = all(withinTol(V_buffered(:,1),V_buffered,pgon.TOL * 21),1);

        % read and remove vertices to be merged
        V_merge = V_buffered(:,idx);
        V_buffered(:,idx) = [];

        % merge vertices to obtain true (unbuffered) vertex
        V_true(:,end+1) = mean(V_merge,2);

        if size(V_true,2) > 2
            % too many points; assume polygon was not buffered
            break
        end
    end

    if isempty(V_buffered)
        % all points were merged, return true vertices
        V = V_true;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
