function pgon = minus(pgon, subtrahend)
% minus - translation of a polygon P by a vector
%
% Syntax:
%    pgon = minus(pgon, subtrahend)
%
% Inputs:
%    pgon - polygon
%    subtrahend - numeric
%
% Outputs:
%    pgon - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       28-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% different types of sets
if isa(pgon, 'polygon') && isnumeric(subtrahend)
    % P - v
    pgon = pgon + -subtrahend;

elseif isnumeric(pgon) && isa(subtrahend, 'polygon')
    % v - P
    pgon = pgon + -subtrahend;

elseif isa(pgon, 'polygon') && isa(subtrahend, 'polygon')
    % P1 - P2
    % throw error
    classname = class(S);
    throw(CORAerror('CORA:notSupported', ...
        sprintf(['The function ''minus'' is not implemented for the class %s except for vectors as a subtrahend.\n', ...
        'If you require to compute the Minkowski difference, use ''minkDiff'' instead.'], classname)));

else
    throw(CORAerror('CORA:noops', pgon, subtrahend));
end

end

% ------------------------------ END OF CODE ------------------------------
