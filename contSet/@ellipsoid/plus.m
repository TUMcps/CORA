function E = plus(E,S,varargin)
% plus - Overloaded '+' operator for approximating the Minkowski sum of an
%    ellipsoid and another set 
%
% Syntax:
%    E = plus(E,S)
%    E = plus(E,S,mode)
%    E = plus(E,S,mode,L)
%
% Inputs:
%    E - ellipsoid object
%    S - set representation/double matrix
%    mode - (optional) type of approximation
%               'inner'
%               'outer':
%               'outer:halder': available when L is empty
%    L - (optional) directions to use for approximation
%
% Outputs:
%    E - ellipsoid object after Minkowski sum
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
% See also: -

% Authors:       Victor Gassmann
% Written:       09-March-2021
% Last update:   04-July-2022 (VG, class array instead of cell array)
%                17-March-2023 (MW, simplify argument pre-processing)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
if nargin > 4
    throw(CORAerror('CORA:tooManyInputArgs',4));
end

% make sure first argument is class argument 
[E,S] = findClassArg(E,S,'ellipsoid');

if isa(S,'ellipsoid')
    [mode,L] = setDefaultValues({'outer:halder',zeros(dim(E),0)},varargin);
else
    [mode,L] = setDefaultValues({'outer',zeros(dim(E),0)},varargin);
end

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {S,'att',{'contSet','numeric'}}; ...
                {mode,'str',{'outer','outer:halder','inner'}}; ...
                {L,'att','numeric'}});

% ind mask for which ellipsoids are empty
ind_empty = representsa_(E,'emptySet',eps);
E(ind_empty) = [];
% Minkowski addition with empty set
if representsa_(E,'emptySet',eps)
    return;
elseif representsa_(S,'origin',eps)
    % adding the origin does not change the set...
    return;
elseif representsa_(S,'emptySet',eps)
    E = ellipsoid.empty(dim(E)); return
end

% dimension check
equalDimCheck(E,S);

% check arguments
if strcmp(mode,'outer:halder') && ~isa(S,'ellipsoid')
    throw(CORAerror('CORA:notSupported',['The method ''outer:halder'' ',...
        'is only implemented for the Minkowski sum of two ellipsoids.']))
end

% check for dimension mismatch...
if size(L,1) ~= dim(E)
    throw(CORAerror('CORA:dimensionMismatch',E,L));
end

N = length(S);

%% different Minkowski additions
if isa(S,'double')
    s = sum(S,2);
    E = ellipsoid(E.Q,E.q+s);
    return;
end

if isa(S,'conPolyZono')
    E = S(1) + E; 
    for i=2:N
        E = S(i) + E; 
    end
    return; 
end

if isa(S,'ellipsoid')
   E = plusEllipsoid([E;S(:)],L,mode);
   return;
end

if isa(S,'interval')
    % convert to ellipsoid
    E = E + ellipsoid(S);
    return;
end

if isa(S,'zonotope')
    % convert to interval (then to ellipsoid)
    E = E + interval(S);
    return;
end

% throw error for all other combinations
throw(CORAerror('CORA:noops',E,S));

% ------------------------------ END OF CODE ------------------------------
