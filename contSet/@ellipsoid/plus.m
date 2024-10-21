function S_out = plus(E,S,varargin)
% plus - Overloaded '+' operator for approximating the Minkowski sum of an
%    ellipsoid and another set or point
%
% Syntax:
%    S_out = E + S
%    S_out = plus(E,S)
%    S_out = plus(E,S,mode)
%    S_out = plus(E,S,mode,L)
%
% Inputs:
%    E - ellipsoid object, numeric
%    S - contSet object (or cell array), numeric
%    mode - (optional) type of approximation
%               'inner'
%               'outer':
%               'outer:halder': available when L is empty
%    L - (optional) directions to use for approximation
%
% Outputs:
%    S_out - set after Minkowski sum
%
% Example: 
%    E1 = ellipsoid(eye(2),[1;-1]);
%    E2 = ellipsoid(diag([1,2]));
%    Ep = E1 + E2;
%    figure; hold on
%    plot(E1); plot(E2);
%    plot(Ep,[1,2],'r');
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal
%       toolbox (ET). In Proceedings of the 45th IEEE Conference on
%       Decision and Control (pp. 1498-1503). IEEE.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_plusEllipsoid

% Authors:       Victor Gassmann
% Written:       09-March-2021
% Last update:   04-July-2022 (VG, class array instead of cell array)
%                17-March-2023 (MW, simplify argument pre-processing)
%                05-October-2024 (MW, remove class array)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,4);

% ensure that numeric is second input argument
[E,S] = reorderNumeric(E,S);

% default values
if isa(S,'ellipsoid')
    [mode,L] = setDefaultValues({'outer:halder',zeros(dim(E),0)},varargin);
else
    [mode,L] = setDefaultValues({'outer',zeros(dim(E),0)},varargin);
end

% check input arguments
inputArgsCheck({{E,'att','ellipsoid'};
                {S,'att',{'cell','contSet','numeric'}}; ...
                {mode,'str',{'outer','outer:halder','inner'}}; ...
                {L,'att','numeric'}});

% Minkowski addition with empty set
if representsa_(E,'emptySet',eps) || (~iscell(S) && representsa_(S,'emptySet',eps))
    S_out = ellipsoid.empty(dim(E));
    return
elseif ~iscell(S) && representsa_(S,'origin',eps)
    % adding the origin does not change the set...
    S_out = E;
    return
end

% dimension checks
equalDimCheck(E,S);
equalDimCheck(E,L);

% addition of vector
if isnumeric(S) && iscolumn(S)
    S_out = ellipsoid(E.Q, E.q+S);
    return;
end

if isa(S,'ellipsoid')
    S_out = priv_plusEllipsoid({E,S},L,mode);
    return;
end

if isa(S,'interval')
    % convert to ellipsoid
    S_out = priv_plusEllipsoid({E,ellipsoid(S)},L,'outer:halder');
    return;
end

if isa(S,'zonotope')
    % convert to ellipsoid via interval enclosure
    S_out = priv_plusEllipsoid({E,ellipsoid(interval(S))},L,'outer:halder');
    return;
end

if isa(S,'conPolyZono')
    S_out = S + E;
    return; 
end

% all supported Minkowski sums convert second input argument to an
% ellipsoid (see above)
if iscell(S) && all(cellfun(@(S_i) isa(S_i,'ellipsoid') || isa(S_i,'interval') || isa(S_i,'zonotope'),'UniformOutput',true))
    S = cellfun(@(S_i) ellipsoid(S_i),S,'UniformOutput',false);
    S_out = priv_plusEllipsoid([E;S],L,mode);
    return
end

% throw error for all other combinations
throw(CORAerror('CORA:noops',E,S));

% ------------------------------ END OF CODE ------------------------------
