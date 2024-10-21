function S_out = or(E,S,varargin)
% or - overloads '|' operator to compute an inner/outer approximation of
%    the union of an ellipsoid and another set representation
%
% Syntax:
%    S_out = E | S
%    S_out = or(E,S)
%    S_out = or(E,S,mode)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object, numeric, cell-array
%    mode - type of approximation (optional):
%               'inner' (inner-approximation)
%               'outer' (outer-approximation)
%
% Outputs:
%    S_out - ellipsoid object
%
% Example: 
%    E1 = ellipsoid([3 -1; -1 1],[1;0]);
%    E2 = ellipsoid([5 1; 1 2],[1;-1]);
%    E = E1 | E2;
%
% References:
%    S. Boyd and L. Vandenberghe, "Convex Optimization", 2004.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       09-March-2021 
% Last update:   15-March-2021
%                04-July-2022 (VG, class array; fixed small bug)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
mode = setDefaultValues({'outer'},varargin);

% ensure that numeric is second input argument
[E,S] = reorderNumeric(E,S);

% check dimensions
equalDimCheck(E,S);

% check input arguments
inputArgsCheck({{E,'att','ellipsoid'};
                {S,'att',{'contSet','numeric','cell'}};
                {mode,'str',{'outer','inner'}}});

% call function with lower precedence
if isa(S,'contSet') && S.precedence < E.precedence
    S_out = or(S,E,varargin{:});
    return
end

% empty set case
if ~representsa_(E,'emptySet',eps) && (~iscell(S) && representsa_(S,'emptySet',eps))
    S_out = E;
    return;
end

% trivial inner approximation: return E
if strcmp(mode,'inner')
    S_out = E;
    return
end

% all to one cell-array
if ~iscell(S)
    S = {S};
end
if any(cellfun(@(S_i) ~isa(S_i,'ellipsoid') && ~isa(S_i,'polytope') && ...
    ~(isnumeric(S_i) && iscolumn(S_i)), S, 'UniformOutput',true))
    throw(CORAerror('CORA:noops',E,S,mode));
end
E_cell = [{E}, cellfun(@(S_i) aux_convert(S_i),S,'UniformOutput',false)];
S_out = priv_orEllipsoidOA(E_cell);

end


% Auxiliary functions -----------------------------------------------------

function E = aux_convert(S)
% helper function to convert multiple operands to ellipsoids correctly
if isnumeric(S) && iscolumn(S)
    E = ellipsoid.origin(length(S)) + S;
elseif isa(S,'polytope')
    E = ellipsoid(S,'outer');
else
    E = S;
end

end

% ------------------------------ END OF CODE ------------------------------
