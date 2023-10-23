function S = minus(S,p)
% minus - translation of a set by a vector
%
% Syntax:
%    S = minus(S,p)
%
% Inputs:
%    S - contSet object
%    p - vector
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: minkDiff

% Authors:       Mark Wetzlinger
% Written:       02-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isnumeric(p)
    % subtrahend is numeric, call 'plus' with negated vector
    S = S + (-p);

elseif isnumeric(S) && isa(p,'contSet')
    % minuend is a vector, subtrahend is a set
    S = -p + S;

else
    % throw error
    classname = class(S);
    throw(CORAerror('CORA:notSupported',...
        sprintf(['The function ''minus'' is not implemented for the class %s except for vectors as a subtrahend.\n', ...
        'If you require to compute the Minkowski difference, use ''minkDiff'' instead.'],classname)));
end

% ------------------------------ END OF CODE ------------------------------
