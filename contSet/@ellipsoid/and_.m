function E = and_(E,S,mode)
% and_ - overloads '&' operator to compute the intersection an ellipsoid
%    and another set representation
%
% Syntax:
%    E = and_(E,S)
%    E = and_(E,S,mode)
%
% Inputs:
%    E              - ellipsoid object
%    S              - set representation (array)
%    mode(optional) - approximation mode ('inner','outer')
%
% Outputs:
%    E - ellipsoid object
%
% Example:
%    E1 = ellipsoid([3 -1; -1 1],[1;0]);
%    E2 = ellipsoid([5 1; 1 2],[1;-1]);
%    E3 = ellipsoid([0.6 -0.4; -0.4 2.2],[0.5;0]);
%    Eo = and(E1,[E2,E3],'outer');
%    Ei = and(E1,[E2,E3],'inner');
%    figure; hold on;
%    plot(E1); plot(E2); plot(E3);
%    plot(Eo); plot(Ei);
%
% References: 
%   [1] A. Kurzhanski et al. "Ellipsoidal Toolbox Manual", 2006
%       https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%            
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   15-October-2019
%                15-March-2021
%                04-July-2022 (VG, replaced cell arrays by class arrays)
% Last revision: 27-March-2023 (MW, rename and_)

% ------------------------------ BEGIN CODE -------------------------------

N = length(S);
% if only center remains
if rank(E)==0 && isa(S,'contSet')
    % if double, already taken care if in 93
    if ~ismethod(S,'contains')
        throw(CORAerror('CORA:noops',E,S{1}));
    end
    if all(arrayfun(@(ii)contains_(S(ii),E.q,'exact',0),1:N))
        E = ellipsoid(zeros(dim(E)),E.q);
    else
        E = ellipsoid;
    end
    return;
end

%% different intersections

% ellipsoid and point
if isa(S,'double')
    % if not all points are equal, overall intersection is empty
    if ~all(all(withinTol(S,repmat(S(:,1),1,size(S,2)),E.TOL))) || ...
        ~contains_(E,S(:,1),'exact',0)
        E = ellipsoid;
    else
        E = ellipsoid(zeros(size(E.Q)),S(:,1));
    end
    return;
end

% ellipsoid and conPolyZono
if isa(S,'conPolyZono')
    if strcmp(mode,'outer')
        E = and_(S(1),E,'exact');
        for i=2:N
            E = and_(S(i),E,'exact');
        end
    else
        throw(CORAerror('CORA:noops',E,S));
    end
    return;
end

% ellipsoid and ellipsoid
if isa(S,'ellipsoid')
    if strcmp(mode,'outer')
        E = andEllipsoidOA(E,S(1));
        for i=2:N
            if representsa_(E,'emptySet',eps)
                break;
            end
            E = andEllipsoidOA(E,S(i));
        end
    else
        E = andEllipsoidIA(E,S);
    end
    return;
end

% ellipsoid and conHyperplane
if isa(S,'conHyperplane')
    for i=1:N
        if representsa_(S(i),'hyperplane',eps)
            E = andHyperplane(E,S(i));
        else
            E = and_(E,polytope(S(i)),mode);
        end
    end
    return;
end

% ellipsoid and polytope
if isa(S,'polytope')
    E = andPolytope(E,S(1),mode);
    for i=2:N
        E = andPolytope(E,S(i),mode);
    end
    return;
end

% ellipsoid and halfspace
if isa(S,'halfspace')
    E = andHalfspace(E,S(1),mode);
    for i=2:N
        E = andHalfspace(E,S(i),mode);
    end
    return;
end

% throw error for remaining combinations
throw(CORAerror('CORA:noops',E,S));

% ------------------------------ END OF CODE ------------------------------
