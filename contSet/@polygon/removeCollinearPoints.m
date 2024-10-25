function pgon = removeCollinearPoints(pgon, varargin)
% removeCollinearPoints - removes collinear points up to a user-defined 
%     tolerance by checking the rank of two vectors computed from two 
%     subsequent pairs of vertices in the list
%
% Syntax:
%    pgon = removeCollinearPoints(pgon,varargin)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    pgon - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       27-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default tolerance
tol = setDefaultValues({0}, varargin);

% read out ordered vertices
V = vertices_(pgon);

% remove collinear vertices
V = removeCollinearVertices2D(V, tol);

% init polygon
pgon = polygon(V);

end

% ------------------------------ END OF CODE ------------------------------
