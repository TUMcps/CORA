function pgon = plus(pgon, summand)
% plus - compute the minkowski sum
%
% Syntax:
%    pgon = plus(pgon, summand)
%
% Inputs:
%    pgon - polygon
%    summand - polygon or numeric
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

% get polygon object
[pgon, summand] = findClassArg(pgon, summand, 'polygon');

% different types of sets
if isnumeric(summand)
    if isscalar(summand)
        % change to vector
        summand = [summand; summand];
    end

    % translate the polygon
    pgon.set = translate(pgon.set, summand');

elseif isa(summand, 'polygon')
    % compute Minkowski sum
    w = warning();
    warning('off');

    T1 = triangulation(pgon);
    T2 = triangulation(summand);
    pgon = [];

    for i = 1:length(T1)
        for j = 1:length(T2)
            V1 = vertices_(T1{i});
            V2 = vertices_(T2{j});

            % add each vertex of V1 to V2
            V1 = reshape(V1,2,1,[]);
            V = V1 + V2;
            V = reshape(V,2,[]);

            % combine
            pgon = pgon | convHull(polygon(V));
        end
    end

    warning(w);
else
    throw(CORAerror('CORA:noops', pgon, summand));
end

end

% ------------------------------ END OF CODE ------------------------------
