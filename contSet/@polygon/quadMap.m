function pgon_2 = quadMap(varargin)
% quadMap - compute tight enclosure of the quadratic map
%
% Syntax:
%    pgon = quadMap(pgon,Q)
%    pgon = quadMap(pgon1,pgon2,Q)
%
% Inputs:
%    pgon,pgon1,pgon2 - polygon
%    Q - cell
%
% Outputs:
%    pgon - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


narginchk(2, 3);

if nargin == 2

    % parse input arguments
    pgon = varargin{1};
    Q = varargin{2};

    if ~all(size(Q) == [2, 1])
        throw(CORAerror('CORA:wrongValue', 'second', "be of the size [2,1]"));
    end

    % compute triangulation of polygon and compute quadratic map
    % for each triangle using polynomial zonotopes
    list = triangulation(pgon);
    pgon_2 = [];

    for i = 1:length(list)
        % convert current triangle to polynomial zonotope
        poly_i = polytope(list{i});
        pZ_i = polyZonotope(poly_i);

        % compute quadratic map
        pZ_i_2 = quadMap(pZ_i, Q);

        % compute polygon enclosure of the result
        pgon_i_2 = polygon(pZ_i_2);

        % unite with results for other triangles
        pgon_2 = pgon_2 | pgon_i_2;
    end

elseif nargin == 3

    % parse input arguments
    pgon1 = varargin{1};
    pgon_i_2 = varargin{2};
    Q = varargin{3};

    if ~all(size(Q) == [2, 1])
        throw(CORAerror('CORA:wrongValue', 'second', "be of the size [2,1]"));
    end

    % compute triangulation of polygons and compute quadratic map
    % for each triangle combinations using polynomial zonotopes
    list1 = triangulation(pgon1);
    list2 = triangulation(pgon_i_2);

    pgon_2 = [];

    for i = 1:length(list1)

        % convert current triangle to polynomial zonotope
        poly_i = polytope(list1{i});
        pZ1 = polyZonotope(poly_i);

        for j = 1:length(list2)

            % convert current triangle to polynomial zonotope
            poly_i = polytope(list2{i});
            pZ_i_2 = polyZonotope(poly_i);

            % compute quadratic map
            pZ_i_2 = quadMap(pZ1, pZ_i_2, Q);

            % compute polygon enclosure of the result
            pZ_i_2 = polygon(pZ_i_2);

            % unite with results for other triangles
            pgon_2 = pgon_2 | pZ_i_2;

        end
    end

end
end

% ------------------------------ END OF CODE ------------------------------
