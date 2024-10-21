function pgon = minkDiff(pgon, subtrahend)
% minkDiff - compute the Minkowski difference
%
% Syntax:
%    pgon = minkDiff(pgon, subtrahend)
%
% Inputs:
%    pgon - polygon, minuend
%    subtrahend - polygon
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

% parse input
narginchk(2, 2);

inputArgsCheck({ ...
    {pgon, 'att', 'polygon'}; ...
    {subtrahend, 'att', {'polygon', 'numeric'}}; ...
    })

% numeric case
if isnumeric(subtrahend)
    pgon = pgon - subtrahend;
    return
end

% compute Mink. diff. by translating the mirrored subtrahend along
% the boundary of the polygon

% shift by center
c = center(subtrahend);
subtrahend = subtrahend - c;
pgon = pgon - c;

% mirror subtrahend
subtrahend = -subtrahend;

% get boundary of pgon
V = vertices_(pgon);

% translate subtrahend along boundary
diff = [];
boundStart = 1;
vs = size(V, 2);
for i = 1:vs
    % boundaries are splitted by nan values

    % take first point on boundary
    Vi = V(:, i);
    if all(isnan(Vi))
        % region boundary
        continue
    end

    % get subsequent point on boundary
    if i + 1 <= vs
        % take next value
        Vi1 = V(:, i+1);
    else
        % take start point of current boundary
        Vi1 = V(:, boundStart);
    end

    % check if end of current boundary is reached
    if all(isnan(Vi1))
        % take start point of current boundary
        Vi1 = V(:, boundStart);
        % shift to start of next boundary
        boundStart = i + 2;
    end

    % linear combination between the two points
    diff = linComb(subtrahend+Vi, subtrahend+Vi1) | diff;
end

pgon.set = subtract(pgon.set, diff.set);

end

% ------------------------------ END OF CODE ------------------------------
