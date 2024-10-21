function P = polytope(SpS,varargin)
% polytope - conversion to polytope objects (only for 2d/3d)
%
% Syntax:
%    P = polytope(SpS)
%    P = polytope(SpS,approxType,splits)
%
% Inputs:
%    SpS - spectraShadow object
%    approxType - 'outer', 'inner'
%    splits - number of splits for refinement
%
% Outputs:
%    P - polytope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg, Tobias Ladner
% Written:       02-August-2023
% Last update:   15-October-2024 (TL, extracted from spectraShadow/plot + simplified)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,3);
[approxType,splits] = setDefaultValues({'inner',100},varargin);
inputArgsCheck({{SpS,'att','spectraShadow'};
    {approxType,'str',{'outer','inner'}}; ...
    {splits,'att','numeric',{'scalar','integer','nonnegative'}} ...
});
if ~any(dim(SpS) == [1,2,3])
    throw(CORAerror('CORA:notSupported','Polytope conversion is only implemented for 1-, 2-, and 3-dimensional spectraShadow.'))
end

% check if set is empty
if representsa(SpS,'emptySet')
    P = polytope.empty(dim(SpS));
    return
end

% switch approximation type
switch approxType
    case 'outer'
        P = aux_polytopeOuter(SpS,splits);
    case 'inner'
        P = aux_polytopeInner(SpS,splits);
    otherwise
        throw(CORAerror('CORA:wrongValue','third',{'inner','outer'}))
end

end


% Auxiliary functions -----------------------------------------------------

function P = aux_polytopeOuter(SpS,splits)
    % Step by step, we construct the H-polytope encircling the
    % spectrahedral shadow

    % shift by center
    c = center(SpS);
    SpS = SpS - c;

    % check dimensions
    n = dim(SpS);
    dirs = aux_getDirections(n,splits);

    % compute support functions 
    P_A = [];
    P_b = [];
    for i=1:size(dirs,2)
        dir = dirs(:,i);
        [val,~] = SpS.supportFunc(dir);
        if val ~= Inf && val~= -Inf
            P_A = [P_A;dir'];
            P_b = [P_b;val];
        end
    end

    % init polytope
    if isempty(P_A)
        % init full-dimensional set
        P = polytope.Inf(n);
    else
        % init polytope with constraints
        P = polytope(P_A,P_b);

        % shift back by center
        P = P + c;
    end

end

function P = aux_polytopeInner(SpS,splits)
    % Step by step, we construct the V-polytope inside the
    % spectrahedral shadow

    % shift by center
    c = center(SpS);
    SpS = SpS - c;

    % check dimensions
    n = dim(SpS);
    dirs = aux_getDirections(n,splits);

    % compute support function in each direction
    V = [];
    for i=1:size(dirs,2)
        dir = dirs(:,i);
        [~,v] = SpS.supportFunc(dir);
        if ~isempty(v)
            V = [V v];
        end
    end

    % init polytope
    if isempty(V)
        % init full-dimensional set
        P = polytope.Inf(n);
    else
        % init polytope with vertices
        if n > 1
            % compute convHull
            Vt = V';
            k = convhulln(Vt);
            k = unique(k(:,1),'stable');
            V = V(:,k);
        end
        P = polytope(V);

        % shift back by center
        P = P + c;
    end

end

function dirs = aux_getDirections(n,splits)
    % check dimensions
    if n == 1
        dirs = [1 -1];

    elseif n == 2 % 2-dimensional case
        % We uniformly distribute points on the unit circle
        phi = linspace(0, 2*pi, splits);
        dirs = [cos(phi);sin(phi)];

    elseif n == 3 % 3-dimensional case
        % We uniformly distribute points on the unit sphere
        dirs = eq_point_set(3-1,splits);
        
    else
        throw(CORAerror('CORA:notSupported','Polytope conversion is only implemented for 1-, 2-, and 3-dimensional spectraShadow'))
    end
end

% ------------------------------ END OF CODE ------------------------------
