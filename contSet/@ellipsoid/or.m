function E = or(E,S,varargin)
% or - overloads '|' operator to compute an inner-/outer-approximation of
%    the union of an ellipsoid and another set representation
%
% Syntax:
%    E = or(E,S)
%    E = or(E,S,mode)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object (or array)
%    mode - type of approximation (optional):
%               'inner' (inner-approximation)
%               'outer' (outer-approximation)
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    E1 = ellipsoid([3 -1; -1 1],[1;0]);
%    E2 = ellipsoid([5 1; 1 2],[1;-1]);
%    E = E1 | E2;
%
% References:
%    Convex Optimization; Boyd
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

%% parsing and checking
% make sure first argument is class argument 
[E,S] = findClassArg(E,S,'ellipsoid');
mode = setDefaultValues({'outer'},varargin);

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {S,'att',{'contSet','numeric','cell'}};
                {mode,'str',{'outer','inner'}}});

% empty set case
if ~representsa_(E,'emptySet',eps) && all(representsa_(S,'emptySet',eps))
    return;
end

% check if equal dims
equalDimCheck(E,S);


if all(representsa_(S,'emptySet',eps))
    return;
end

%% different unions
if isa(S,'double')
    V = S;
    if strcmp(mode,'outer')
        E_p = ellipsoid.enclosePoints(V,'min-vol');
        E = orEllipsoidOA([E,E_p]);
        return;
    else
        P = polytope.enclosePoints(V);
        E = or(E,P,'inner');
        return;
    end
end

if isa(S,'ellipsoid')
    if strcmp(mode,'outer')
        E = orEllipsoidOA([E;S(:)]);
        return;
    else
        throw(CORAerror('CORA:noops',E,S,'inner'));
    end
end

if isa(S,'polytope')
    if strcmp(mode,'outer')
        E_S = ellipsoid.array(1,numel(S));
        for i=1:numel(S)
            E_S(i) = ellipsoid(S,'outer');
        end
        E = orEllipsoid([E,E_S]);
    else
        % not implemented
        %E = orEllipsoidIA([E,ellipsoid(S,'inner')]);
        throw(CORAerror('CORA:noops',E,S));
    end
end

% throw error for all other combinations
throw(CORAerror('CORA:noops',E,S));

% ------------------------------ END OF CODE ------------------------------
